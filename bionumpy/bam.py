from itertools import accumulate, repeat, takewhile
from .npdataclassstream import streamable
from .datatypes import Interval, Bed6, BamEntry
from .file_buffers import FileBuffer
from .cigar import count_reference_length, split_cigar, CigarOpArray
from .datatypes import SAMEntry
from npstructures.raggedshape import RaggedView
from npstructures import RaggedArray, npdataclass
import numpy as np
from .sequences import as_sequence_array, Sequence, create_sequence_array_from_already_encoded_data
from .encodings import AlphabetEncoding, BaseEncoding, QualityEncoding, Encoding, NumericEncoding
from .encodings.alphabet_encoding import BamArray, CigarOpArray
# from bionumpy.bnpdataclass import bnpdataclass

# BamEncoding = AlphabetEncoding("=ACMGRSVTWYHKDBN")
# BamArray = get_alphabet_array_class("=ACMGRSVTWYHKDBN")


class BamBuffer(FileBuffer):
    dataclass = BamEntry

    def __init__(self, data, delimiters, header_data):
        super().__init__(data, delimiters)
        self._chromosome_names = as_sequence_array([header[0] for header in header_data])
        self._data = np.asarray(self._data)

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
        new_start = lambda start, _: start + int.from_bytes(chunk[start:start+4], byteorder="little") + 4
        _starts = accumulate(repeat(None), new_start, initial=0)
        starts = list(takewhile(lambda start: start <= len(chunk), _starts))
        return starts

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
        ref_id = self._get_ints(4, 4, np.int32)
        pos = self._get_ints(8, 4, np.int32)
        chromosome = self._chromosome_names[ref_id]
        l_read_name = self._data[self._new_lines+12]
        mapq = self._data[self._new_lines+13]
        l_seq = self._get_ints(20, 4, np.int32)
        n_seq_bytes = (l_seq+1)//2
        n_cigar_op = self._get_ints(16, 2, np.uint16)
        flag = self._get_ints(18, 1, np.uint8)
        n_cigar_bytes = n_cigar_op*4
        read_names = self._move_intervals_to_ragged_array(self._new_lines+36, self._new_lines+36+l_read_name-1)
        cigars = self._move_intervals_to_ragged_array(self._new_lines+36+l_read_name,
                                                      self._new_lines+36+l_read_name+n_cigar_bytes, as_sequence=False)

        cigars = RaggedArray(cigars.ravel().view(np.uint32), cigars.shape.lengths//4)
        cigar_cymbol, cigar_length = split_cigar(cigars)
        sequences = self._move_intervals_to_ragged_array(self._new_lines+36+l_read_name+n_cigar_bytes,
                                                         self._new_lines+36+l_read_name+n_cigar_bytes+n_seq_bytes)
        sequences = (((sequences.ravel()[:, None]) >> (4*np.arange(2, dtype=np.uint8))).ravel() & np.uint8(15)).view(BamArray)
        new_sequences = RaggedArray(sequences, n_seq_bytes*2)
        view = RaggedView(new_sequences.shape.starts, l_seq)
        new_sequences = new_sequences[view]
        quals = self._move_intervals_to_ragged_array(self._new_lines+36+l_read_name+n_cigar_bytes+n_seq_bytes,
                                                     self._new_lines+36+l_read_name+n_cigar_bytes+n_seq_bytes+l_seq)# +33
        return BamEntry(chromosome, read_names, flag, pos, mapq, cigar_cymbol, cigar_length, new_sequences, quals)

    def count_entries(self):
        return len(self._new_lines)


class BamIntervalBuffer(BamBuffer):
    dataclass = Bed6

    def get_data(self):
        ref_id = self._get_ints(4, 4, np.int32)
        pos = self._get_ints(8, 4, np.int32)
        chromosome = self._chromosome_names[ref_id]
        l_read_name = self._data[self._new_lines+12]
        mapq = self._data[self._new_lines+13]
        n_cigar_op = self._get_ints(16, 2, np.uint16)
        flag = self._get_ints(18, 1, np.uint8)
        n_cigar_bytes = n_cigar_op*4
        read_names = self._move_intervals_to_ragged_array(self._new_lines+36, self._new_lines+36+l_read_name-1)
        cigars = self._move_intervals_to_ragged_array(self._new_lines+36+l_read_name,
                                                      self._new_lines+36+l_read_name+n_cigar_bytes, as_sequence=False)

        cigars = RaggedArray(cigars.ravel().view(np.uint32), cigars.shape.lengths//4)
        cigar_cymbol, cigar_length = split_cigar(cigars)

        strand = flag & np.uint16(16)
        strand = np.where(strand, ord("-"), ord("+"))[:, None].view(Sequence)
        strand.encoding = BaseEncoding
        length = count_reference_length(cigar_cymbol, cigar_length)
        return Bed6(chromosome,
                    pos,
                    pos+length,
                    read_names,
                    mapq,
                    strand)


@streamable()
def alignment_to_interval(alignment):
    strand = alignment.flag & np.uint16(16)
    strand = np.where(strand, ord("-"), ord("+"))[:, None].view(Sequence)
    strand.encoding = BaseEncoding
    length = count_reference_length(alignment.cigar_op, alignment.cigar_length)
    return Bed6(alignment.chromosome,
                alignment.position,
                alignment.position+length,
                alignment.name,
                alignment.mapq,
                strand)
