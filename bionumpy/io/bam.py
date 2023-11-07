from functools import lru_cache
from itertools import accumulate, repeat, takewhile, chain
import numpy as np
from npstructures.raggedshape import RaggedView
from npstructures import RaggedArray, ragged_slice

from ..datatypes import Bed6, BamEntry
from ..alignments.cigar import count_reference_length, split_cigar
from ..encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from ..encodings import BaseEncoding
from ..encodings.alphabet_encoding import BamEncoding
from .file_buffers import FileBuffer
from ..util import cached_property


class BamBufferExtractor:
    def __init__(self, data, new_entries, chromosome_names):
        self._data = data
        self._new_lines = new_entries
        self._chromosome_names = chromosome_names
        self._functions = [self._get_chromosome,
                           self._get_read_name,
                           self._get_flag,
                           self._get_position,
                           self._get_mapq,
                           self._get_cigar_symbol,
                           self._get_cigar_length,
                           self._get_sequences,
                           self._get_quality]

    def _get_ints(self, offsets, n_bytes, dtype):
        tmp = self._data[(self._new_lines+offsets)[:, None] + np.arange(n_bytes)].ravel()
        ints = (tmp).view(dtype).ravel()
        assert len(ints) == len(self._new_lines), (len(ints), offsets, len(self._new_lines), n_bytes, dtype)
        return ints

    def _get_cigar_bytes(self):
        n_cigar_op = self._get_ints(16, 2, np.uint16)
        n_cigar_bytes = n_cigar_op * 4
        return n_cigar_bytes

    def _get_quality(self):
        return ragged_slice(self._data, self._quality_start, self._quality_start + self._get_sequence_length())

    def _get_sequences(self):
        l_seq = self._get_sequence_length()
        n_seq_bytes = (l_seq + 1) // 2
        sequences = ragged_slice(self._data, self._sequence_start, self._quality_start)
        sequences = EncodedArray(
            (((sequences.ravel()[:, None]) >> (4 * np.arange(2, dtype=np.uint8)[::-1])).ravel() & np.uint8(15)),
            BamEncoding)
        new_sequences = EncodedRaggedArray(sequences, n_seq_bytes * 2)
        view = RaggedView(new_sequences._shape.starts, l_seq)
        new_sequences = new_sequences[view]
        return new_sequences

    def _get_cigar(self):
        cigars = ragged_slice(self._data, self._cigar_start, self._sequence_start)
        cigars = RaggedArray(cigars.ravel().view(np.uint32), cigars.lengths // 4)
        cigar_symbol, cigar_length = split_cigar(cigars)
        return cigar_symbol, cigar_length

    def _get_cigar_symbol(self):
        return self._get_cigar()[0]

    def _get_cigar_length(self):
        return self._get_cigar()[1]

    @cached_property
    def _read_name_start(self):
        return self._new_lines + 36

    @cached_property
    def _cigar_start(self):
        return self._read_name_start + self._get_read_name_length()

    @cached_property
    def _sequence_start(self):
        return self._cigar_start + self._get_cigar_bytes()

    @cached_property
    def _quality_start(self):
        return self._sequence_start + (self._get_sequence_length() + 1) // 2

    def _get_read_name(self):
        read_names = ragged_slice(self._data, self._read_name_start, self._cigar_start-1)
        read_names = EncodedRaggedArray(
            EncodedArray(read_names.ravel(), BaseEncoding), read_names.shape)
        return read_names

    def _get_flag(self):
        return self._get_ints(18, 2, np.uint16)

    @lru_cache(None)
    def _get_sequence_length(self):
        return self._get_ints(20, 4, np.int32)

    def _get_read_name_length(self):
        return self._data[self._new_lines + 12]

    def _get_mapq(self):
        return self._data[self._new_lines + 13]

    def _get_position(self):
        return self._get_ints(8, 4, np.int32)

    def _get_chromosome(self):
        ref_id = self._get_ints(4, 4, np.int32)
        chromosome = self._chromosome_names[ref_id]
        return chromosome

    def get_field_by_number(self, i):
        return self._functions[i]()



class BamBuffer(FileBuffer):
    '''
    https://samtools.github.io/hts-specs/SAMv1.pdf
    '''
    dataclass = BamEntry

    def __init__(self, data, delimiters, header_data):
        delimiters = np.asanyarray(delimiters)
        super().__init__(data, delimiters)
        self._chromosome_names = as_encoded_array([header[0] for header in header_data])
        self._buffer_extractor = BamBufferExtractor(data, delimiters, self._chromosome_names)
        self._data = np.asarray(self._data)

    @property
    def data(self):
        return self._data

    @property
    def n_lines(self):
        return len(self._new_lines)

    @classmethod
    def _read_int(self, file_object):
        return int.from_bytes(file_object.read(4), byteorder="little")

    @classmethod
    def _read_zero_term(cls, file_object):
        chars = []
        while True:
            chars.append(file_object.read(1))
            if chars[-1] == b"\x00":
                break
        return "".join(c.decode("ascii") for c in chars[:-1])

    @classmethod
    def read_header(cls, file_object):
        magic = file_object.read(4)
        assert magic == b"BAM\1", magic
        header_length = cls._read_int(file_object)
        file_object.read(header_length)
        n_ref = cls._read_int(file_object)
        return cls._handle_refs(n_ref, file_object)

    @classmethod
    def _handle_refs(cls, n_refs, file_object):
        info = []
        for _ in range(n_refs):
            ref_n = cls._read_int(file_object)
            name = cls._read_zero_term(file_object)
            sequence_length = cls._read_int(file_object)
            info.append((name, sequence_length))
        return info

    @staticmethod
    def _find_starts(chunk):
        chunk = bytes(chunk)
        def new_start(start, _):
            return start + int.from_bytes(chunk[start:start+4], byteorder="little") + 4
        # new_start = lambda start, _: start + int.from_bytes(chunk[start:start+4], byteorder="little") + 4
        _starts = accumulate(repeat(0), new_start)# chain([0], accumulate(repeat(0), new_start))
        starts = list(takewhile(lambda start: start <= len(chunk), _starts))
        return starts

    @classmethod
    def contains_complete_entry(cls, chunks):
        return True

    @classmethod
    def from_raw_buffer(cls, chunk, header_data):
        chunk = np.asarray(chunk)
        starts = cls._find_starts(chunk)
        return cls(chunk[:starts[-1]], starts[:-1], header_data)

    def _get_ints(self, offsets, n_bytes, dtype):
        tmp = self._data[(self._new_lines+offsets)[:, None] + np.arange(n_bytes)].ravel()
        ints = (tmp).view(dtype).ravel()
        assert len(ints) == len(self._new_lines), (len(ints), offsets, len(self._new_lines), n_bytes, dtype)
        return ints

    def get_data(self):
        l_read_name = self._get_read_name_length()
        l_seq = self._get_sequence_length()
        n_seq_bytes = (l_seq+1)//2
        n_cigar_bytes = self._get_cigar_bytes()
        read_name_start = self._new_lines + 36
        cigar_start = read_name_start + l_read_name
        read_name_end = cigar_start - 1
        sequence_start = cigar_start + n_cigar_bytes
        cigar_cymbol, cigar_length = self._get_cigar(cigar_start, sequence_start)
        quality_start = sequence_start + n_seq_bytes
        return BamEntry(self._buffer_extractor.get_field_by_number(0),
                        self._buffer_extractor.get_field_by_number(1),
                        self._buffer_extractor.get_field_by_number(2),
                        self._buffer_extractor.get_field_by_number(3),
                        self._buffer_extractor.get_field_by_number(4),
                        self._buffer_extractor.get_field_by_number(5),
                        self._buffer_extractor.get_field_by_number(6),
                        self._buffer_extractor.get_field_by_number(7),
                        self._buffer_extractor.get_field_by_number(8))
                        # self._get_sequences(l_seq, n_seq_bytes, quality_start, sequence_start),
                        # self._get_quality(l_seq, quality_start))

    def _get_cigar_bytes(self):
        n_cigar_op = self._get_ints(16, 2, np.uint16)
        n_cigar_bytes = n_cigar_op * 4
        return n_cigar_bytes

    def count_entries(self):
        return len(self._new_lines)

    def _get_quality(self, l_seq, quality_start):
        return ragged_slice(self._data, quality_start, quality_start + l_seq)

    def _get_sequences(self, l_seq, n_seq_bytes, quality_start, sequence_start):
        sequences = ragged_slice(self._data, sequence_start, quality_start)
        sequences = EncodedArray(
            (((sequences.ravel()[:, None]) >> (4 * np.arange(2, dtype=np.uint8)[::-1])).ravel() & np.uint8(15)),
            BamEncoding)
        new_sequences = EncodedRaggedArray(sequences, n_seq_bytes * 2)
        view = RaggedView(new_sequences._shape.starts, l_seq)
        new_sequences = new_sequences[view]
        return new_sequences

    def _get_cigar(self, cigar_start, sequence_start):
        cigars = ragged_slice(self._data, cigar_start, sequence_start)
        cigars = RaggedArray(cigars.ravel().view(np.uint32), cigars.lengths // 4)
        cigar_symbol, cigar_length = split_cigar(cigars)
        return cigar_symbol, cigar_length

    def _get_read_name(self, read_name_end, read_name_start):
        read_names = ragged_slice(self._data, read_name_start, read_name_end)
        read_names = EncodedRaggedArray(
            EncodedArray(read_names.ravel(), BaseEncoding), read_names.shape)
        return read_names

    def _get_flag(self):
        return self._get_ints(18, 2, np.uint16)

    def _get_sequence_length(self):
        return self._get_ints(20, 4, np.int32)

    def _get_read_name_length(self):
        return self._data[self._new_lines + 12]

    def _get_mapq(self):
        return self._data[self._new_lines + 13]

    def _get_position(self):
        return self._get_ints(8, 4, np.int32)

    def _get_chromosome(self):
        ref_id = self._get_ints(4, 4, np.int32)
        chromosome = self._chromosome_names[ref_id]
        return chromosome


class BamIntervalBuffer(BamBuffer):
    dataclass = Bed6

    def get_data(self):
        ref_id = self._get_ints(4, 4, np.int32)
        pos = self._get_ints(8, 4, np.int32)
        chromosome = self._chromosome_names[ref_id]
        l_read_name = self._data[self._new_lines+12]
        mapq = self._data[self._new_lines+13]
        n_cigar_op = self._get_ints(16, 2, np.uint16)
        flag = self._get_ints(18, 2, np.uint16)
        n_cigar_bytes = n_cigar_op*4
        read_names = ragged_slice(self._data, self._new_lines+36, self._new_lines+36+l_read_name-1)
        read_names = EncodedRaggedArray(
            EncodedArray(read_names.ravel(), BaseEncoding), read_names.shape)
        cigars = ragged_slice(self._data, self._new_lines+36+l_read_name,
                              self._new_lines+36+l_read_name+n_cigar_bytes)

        cigars = RaggedArray(cigars.ravel().view(np.uint32), cigars.lengths//4)
        cigar_cymbol, cigar_length = split_cigar(cigars)

        strand = flag & np.uint16(16)
        strand = EncodedArray(np.where(strand, ord("-"), ord("+"))[:, None], BaseEncoding)
        strand.encoding = BaseEncoding
        length = count_reference_length(cigar_cymbol, cigar_length)
        return Bed6(chromosome,
                    pos,
                    pos+length,
                    read_names,
                    mapq,
                    strand)
