from .npdataclassstream import streamable
from .datatypes import Interval, Bed6
from .file_buffers import FileBuffer
from .cigar import count_reference_length
from .datatypes import SAMEntry, StrandedInterval
from npstructures.raggedshape import RaggedView
from npstructures import RaggedArray, npdataclass
import numpy as np
from .sequences import Sequences, as_sequence_array, Sequence
from .encodings import AlphabetEncoding, BaseEncoding

BamEncoding = AlphabetEncoding("=ACMGRSVTWYHKDBN")


@npdataclass
class Dummy:
    chromosome: str
    name: str
    flag: int
    position: int
    mapq: int
    cigar: str
    sequence: str
    quality: str


class BamBuffer(FileBuffer):
    dataclass = Dummy
    def __init__(self, data, delimiters, header_data):
        super().__init__(data, delimiters)
        self._chromosome_names = as_sequence_array([header[0] for header in header_data])
        #self._header_data = header_data

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

    @classmethod
    def from_raw_buffer(cls, chunk, header_data):
        starts = [0]
        while (starts[-1]+3) < chunk.size-1:
            line_size = (chunk[starts[-1]:starts[-1]+4]).view(np.int32)[0]
            starts.append(starts[-1]+line_size+4)
        if starts[-1] > chunk.size:
            starts.pop()
        return cls(chunk[:starts[-1]], starts[:-1], header_data)

    def _get_ints(self, offsets, n_bytes, dtype):
        tmp = self._data[(self._new_lines+offsets)[:, None] + np.arange(n_bytes)].ravel()
        ints = (tmp).view(dtype).ravel()
        assert len(ints) == len(self._new_lines), (len(ints), self._new_lines)
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
        sequences = self._move_intervals_to_ragged_array(self._new_lines+36+l_read_name+n_cigar_bytes,
                                                         self._new_lines+36+l_read_name+n_cigar_bytes+n_seq_bytes)
        new_sequences = Sequences(((sequences.ravel()[:, None]) >> np.arange(2, dtype=np.uint8)).ravel() & np.uint8(15), n_seq_bytes*2, encoding=BamEncoding)

        quals = self._move_intervals_to_ragged_array(self._new_lines+36+l_read_name+n_cigar_bytes+n_seq_bytes,
                                                     self._new_lines+36+l_read_name+n_cigar_bytes+n_seq_bytes+l_seq)+33
        return Dummy(chromosome, read_names, flag, pos, mapq, cigars, new_sequences, quals)
        return SAMEntry(read_names, pos, )

class BAMReader:
    def handle_alignment(self):
        block_size = self._read_int()
        ref_id = self._read_int()
        pos = self._read_int()
        l_read_name, mapq = (self._read8() for _ in range(2))
        _bin, n_cigar_op, flag = (self._read16() for _ in range(3))
        l_seq = self._read_int()
        next_refID, next_pos, t_len = (self._read_int() for _ in range(3))
        read_name = self._read_zero_term()
        cigar = self._file.read(n_cigar_op)
        seq = self._file.read((l_seq+1)//2)

    def __init__(self, filename):
        if isinstance(filename, str):
            self._file = open(filename, "rb")
        else:
            self._file = filename
        self._info = self._handle_header()

    def _read_int(self):
        return int.from_bytes(self._file.read(4), byteorder="little")

    def _handle_header(self):
        magic = self._file.read(4)
        assert magic == b"BAM\1", magic
        header_length = self._read_int()
        self._file.read(header_length)
        n_ref = self._read_int()
        return self._handle_refs(n_ref)

    def _handle_refs(self, n_refs):
        info = []
        for _ in range(n_refs):
            ref_n = self._read_int()
            name = self._read_zero_term()

            sequence_length = self._read_int()
            info.append((name, sequence_length))
        return info

    def read(self, n):
        self._reader.read(n)

    def _read8(self):
        return int.from_bytes(self._file.read(1), "little")

    def _read16(self):
        return int.from_bytes(self._file.read(2), "little")

    def _read_zero_term(self):
        chars = []
        while True:
            chars.append(self._file.read(1))
            if chars[-1] == b"\x00":
                break
        return "".join(c.decode("ascii") for c in chars[:-1])
        
    def read_block_lengths(self):
        for _ in range(1000):
            block_size = self._read_int()
            if block_size == 0:
                return
            self._file.read(block_size)
            

@streamable
def alignment_to_interval(alignment):
    strand = alignment.flag & np.uint16(16)
    strand = np.where(strand, ord("-"), ord("+"))[:, None].view(Sequence)
    strand.encoding = BaseEncoding
    length = count_reference_length(alignment.cigar)
    return Bed6(alignment.chromosome,
                alignment.position,
                alignment.position+count_reference_length(alignment.cigar),
                alignment.name,
                alignment.mapq,
                strand)
