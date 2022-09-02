from .file_buffers import FileBuffer
from npstructures.raggedshape import RaggedView
from npstructures import RaggedArray
import numpy as np


class BamBuffer(FileBuffer):
    @classmethod
    def from_raw_buffer(cls, chunk):
        starts = [0]
        while (starts[-1]+3) < chunk.size:
            line_size = (chunk[starts[-1]:starts[-1]+4]).view(np.int32)[0]
            print("#", starts[-1], line_size)
            starts.append(starts[-1]+line_size+4)
        return cls(chunk[:starts[-2]], starts[:-2])

    def get_data(self):
        print(self._new_lines)
        indices = self._new_lines[:, None]+8+np.arange(4)
        pos = self._data[indices].view(np.int32)
        l_read_name = self._data[self._new_lines+12]
        l_seq = self._data[self._new_lines[:, None]+20+np.arange(4)].view(np.int32).ravel()
        n_cigar_op = self._data[self._new_lines[:, None]+16+np.arange(2)].view(np.uint16).ravel()
        seq_starts = self._new_lines+36 + l_read_name + n_cigar_op
        print(n_cigar_op.shape)
        print(seq_starts, l_seq)
        view = RaggedView(seq_starts, l_seq)
        indices, shape = view.get_flat_indices()
        seq = RaggedArray(self._data[indices], shape)
        print(seq)
        

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
            
