import numpy as np
from npstructures import RaggedArray, RaggedView
from .encodings import BaseEncoding, ACTGTwoBitEncoding
from .sequences import Sequences
NEWLINE = 10

class FileBuffer:
    _buffer_divisor = 1
    COMMENT = 0

    def __init__(self, data, new_lines):
        self._data = np.asanyarray(data)
        self._new_lines = np.asanyarray(new_lines)
        self._is_validated = False
        self.size = self._data.size

    def _move_intervals_to_2d_array(self, starts, ends, fill_value=0):
        n_intervals = starts.size
        n_chars = (ends-starts)
        from_indices, _ = RaggedView(starts, n_chars).get_flat_indices()
        max_chars = np.max(n_chars)
        array = np.full(n_intervals*max_chars, fill_value, dtype=np.uint8)
        to_indices, _ = RaggedView(max_chars*np.arange(1,n_intervals+1)-n_chars, n_chars).get_flat_indices()
        array[to_indices] = self._data[from_indices]
        return array.reshape((n_intervals, max_chars))

    def _move_intervals_to_ragged_array(self, starts, ends=None, lens=None):
        if lens is None:
            lens = ends-starts
        indices, shape = RaggedView(starts, lens).get_flat_indices()
        return RaggedArray(self._data[indices], shape)

    @classmethod
    def from_raw_buffer(cls, chunk):
        raise NotImplemented

    def validate_if_not(self):
        if not self._is_validated:
            self._validate()

class OneLineBuffer(FileBuffer):
    n_lines_per_entry = 2
    _buffer_divisor = 32
    @classmethod
    def from_raw_buffer(cls, chunk):
        new_lines = np.flatnonzero(chunk==NEWLINE)
        n_lines = new_lines.size
        assert n_lines >= cls.n_lines_per_entry, "No complete entry in buffer"
        new_lines = new_lines[:n_lines-(n_lines%cls.n_lines_per_entry)]
        return cls(chunk[:new_lines[-1]+1], new_lines)

    def get_sequences(self):
        self.validate_if_not()
        sequence_starts = self._new_lines[::self.n_lines_per_entry]+1
        sequence_lens = self._new_lines[1::self.n_lines_per_entry]-sequence_starts
        indices, shape = RaggedView(sequence_starts, sequence_lens).get_flat_indices()
        m = indices.size
        d = m % self._buffer_divisor
        seq = np.empty(m-d+self._buffer_divisor, dtype=self._data.dtype)
        seq[:m] = self._data[indices]
        return Sequences(seq, shape)
    
    def _validate(self):
        n_lines = self._new_lines.size
        assert n_lines % self.n_lines_per_entry == 0, "Wrong number of lines in buffer"
        header_idxs = self._new_lines[self.n_lines_per_entry-1:-1:self.n_lines_per_entry]+1
        assert np.all(self._data[header_idxs]==self.HEADER)
        self._is_validated = True

class TwoLineFastaBuffer(OneLineBuffer):
    HEADER= 62
    n_lines_per_entry = 2
    _encoding = BaseEncoding

class FastQBuffer(OneLineBuffer):
    HEADER= 64
    n_lines_per_entry = 4
    _encoding = BaseEncoding
