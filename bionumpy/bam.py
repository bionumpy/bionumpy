from .file_buffers import FileBuffer
from .datatypes import SAMEntry
from npstructures.raggedshape import RaggedView
from npstructures import RaggedArray, npdataclass
import numpy as np
from .sequences import Sequences
from .encodings import AlphabetEncoding


BamEncoding = AlphabetEncoding("=ACMGRSVTWYHKDBN")


@npdataclass
class Dummy:
    name: str
    flag: int
    position: int
    mapq: int
    cigar: str
    sequence: str
    quality: str


class BamBuffer(FileBuffer):
    @classmethod
    def from_raw_buffer(cls, chunk):
        starts = [0]
        while (starts[-1]+3) < chunk.size:
            line_size = (chunk[starts[-1]:starts[-1]+4]).view(np.int32)[0]
            print("#", starts[-1], line_size)
            starts.append(starts[-1]+line_size+4)
        return cls(chunk[:starts[-1]], starts[:-1])

    def _get_ints(self, offsets, n_bytes, dtype):
        tmp = self._data[(self._new_lines+offsets)[:, None] + np.arange(n_bytes)].ravel()
        print(tmp, tmp.dtype, tmp.size)
        ints = (tmp).view(dtype).ravel()
        assert len(ints) == len(self._new_lines), (len(ints), self._new_lines)
        return ints

    def get_data(self):
        pos = self._get_ints(8, 4, np.int32)
        print("pos", pos)
        l_read_name = self._data[self._new_lines+12]
        mapq = self._data[self._new_lines+13]
        l_seq = self._get_ints(20, 4, np.int32)
        n_seq_bytes = (l_seq+1)//2
        n_cigar_op = self._get_ints(16, 2, np.uint16)
        print(n_cigar_op)
        flag = self._get_ints(18, 1, np.uint8)
        n_cigar_bytes = n_cigar_op*4
        read_names = self._move_intervals_to_ragged_array(self._new_lines+36, self._new_lines+36+l_read_name-1)
        cigars = self._move_intervals_to_ragged_array(self._new_lines+36+l_read_name,
                                                      self._new_lines+36+l_read_name+n_cigar_bytes)
        sequences = self._move_intervals_to_ragged_array(self._new_lines+36+l_read_name+n_cigar_bytes,
                                                         self._new_lines+36+l_read_name+n_cigar_bytes+n_seq_bytes)
        new_sequences = Sequences(((sequences.ravel()[:, None]) >> np.arange(2, dtype=np.uint8)).ravel() & np.uint8(15), n_seq_bytes*2, encoding=BamEncoding)

# Sequences(
# BamEncoding.decode(((sequences.ravel()[:, None]) >> np.arange(2, dtype=np.uint8)).ravel() & np.uint8(15)), n_seq_bytes*2)
        print(">>>", repr(new_sequences._data))
        print(">>>", str(new_sequences._data))
        quals = self._move_intervals_to_ragged_array(self._new_lines+36+l_read_name+n_cigar_bytes+n_seq_bytes,
                                                     self._new_lines+36+l_read_name+n_cigar_bytes+n_seq_bytes+l_seq)+33
        return Dummy(read_names, flag, pos, mapq, cigars, new_sequences, quals)
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
        print(block_size, ref_id, pos, l_read_name, mapq, read_name)

    def __init__(self, filename):
        if isinstance(filename, str):
            self._file = open(filename, "rb")
        else:
            self._file = filename
        self._info = self._handle_header()
        print(self._info)

    def _read_int(self):
        return int.from_bytes(self._file.read(4), byteorder="little")

    def _handle_header(self):
        magic = self._file.read(4)
        assert magic == b"BAM\1", magic
        header_length = self._read_int()
        print(header_length)
        self._file.read(header_length)
        n_ref = self._read_int()
        return self._handle_refs(n_ref)

    def _handle_refs(self, n_refs):
        info = []
        for _ in range(n_refs):
            ref_n = self._read_int()
            name = self._read_zero_term()

            n_sequences = self._read_int()
            info.append((ref_n, name, n_sequences))
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
            print(block_size)
            self._file.read(block_size)
            
