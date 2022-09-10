import cupy as cp
import numpy as np

from ..file_buffers import OneLineBuffer, TwoLineFastaBuffer
from npstructures import RaggedArray, RaggedShape
from .sequences import CPSequences
from ..encodings import BaseEncoding, QualityEncoding
from ..datatypes import SequenceEntry, SequenceEntryWithQuality

NEWLINE = 10

class CPOneLineBuffer(OneLineBuffer):
    @classmethod
    def from_raw_buffer(cls, chunk) -> "OneLineBuffer":
        """Create a buffer with full entries

        Extract complete entries, i. e. a number of lines that is divisible by lines per entry

        Parameters
        ----------
        chunk : np.ndarray
            Raw buffer with data that might end with incomplete entry

        Returns
        -------
        'OneLineBuffer'
            Buffer with complete entries

        Examples
        --------
        8

        """
        new_lines = cp.flatnonzero(chunk == NEWLINE)
        n_lines = new_lines.size
        assert n_lines >= cls.n_lines_per_entry, "No complete entry in buffer"
        new_lines = new_lines[: n_lines - (n_lines % cls.n_lines_per_entry)]
        return cls(chunk[: new_lines[-1] + 1], new_lines)

    def get_sequences(self) -> CPSequences:
        self.validate_if_not()
        sequence_starts = self._new_lines[:: self.n_lines_per_entry] + 1
        sequence_lens = self._new_lines[1 :: self.n_lines_per_entry] - sequence_starts
        indices, shape = RaggedView(sequence_starts, sequence_lens).get_flat_indices()
        m = indices.size
        d = m % self._buffer_divisor
        seq = cp.empty(m - d + self._buffer_divisor, dtype=self._data.dtype)
        seq[:m] = self._data[indices]
        return CPSequences(seq, shape)

    def get_data(self):
        self.validate_if_not()

        #print(self._new_lines.shape)
        #starts = cp.insert(self._new_lines, 0, -1) # TODO fix for cupy
        starts = cp.concatenate((cp.asanyarray([-1]), self._new_lines))

        lengths = cp.diff(starts)
        self.lines = CPSequences(self._data, RaggedShape(lengths))
        sequences = self.lines[1 :: self.n_lines_per_entry, :-1]

        headers = self.lines[:: self.n_lines_per_entry, 1:-1]
        return SequenceEntry(headers, sequences)

class CPTwoLineFastaBuffer(CPOneLineBuffer):
    HEADER = 62
    n_lines_per_entry = 2
    _encoding = BaseEncoding
    dataclass = SequenceEntry
