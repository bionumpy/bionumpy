from .delimited_buffers import DelimitedBuffer
from ..datatypes import BedGraph
from .strops import split
from .exceptions import FormatException
from ..encoded_array import BaseEncoding, EncodedArray
import numpy as np


class WigBuffer(DelimitedBuffer):
    dataclass = BedGraph
    DELIMITER = '\t'

    def __init__(self, data: EncodedArray, new_lines: np.ndarray, delimiters: np.ndarray = None, header_data=None):
        delimiters_mask = data == self.DELIMITER
        delimiters_mask[new_lines] = True
        delimiters = np.append(np.flatnonzero(delimiters_mask), data.size)
        super().__init__(data, new_lines, delimiters, header_data)
        starts, ends = self._calculate_col_starts_and_ends(data, delimiters)
        n_fields = next(i for i, d in enumerate(ends) if data[d] == '\n') + 1
        self._n_cols = n_fields
        self.__col_starts = starts.reshape(-1, n_fields)
        self.__col_ends = ends.reshape(-1, n_fields)

    def _calculate_col_starts_and_ends(self, data, delimiters):
        comment_mask = (data[delimiters[:-1]] == '\n') & (data[delimiters[:-1]+1] == self.COMMENT)
        comment_mask = np.flatnonzero(comment_mask)
        start_delimiters = np.insert(np.delete(delimiters, comment_mask), 0, -1)[:-1]+1
        end_delimiters = np.delete(delimiters, comment_mask+1)
        return start_delimiters, end_delimiters

    def _col_ends(self, col):
        return self.__col_ends[:, col]

    def _col_starts(self, col):
        return self.__col_starts[:, col]

    def _validate(self):
        self._validated = True
        # chunk = self._data
        # delimiters = self._delimiters[1:]
        # n_delimiters_per_line = (
        #     next(i for i, d in enumerate(delimiters) if chunk[d] == '\n') + 1
        # )
        # self._n_cols = n_delimiters_per_line
        # should_be_new_lines = chunk[delimiters[n_delimiters_per_line-1::n_delimiters_per_line]]
        # if delimiters.size % n_delimiters_per_line != 0 or np.any(should_be_new_lines != "\n"):
        #     offending_line = np.flatnonzero(should_be_new_lines != "\n")[0]
        #     lines = split(self._data, '\n')
        #     raise FormatException(f"Irregular number of delimiters per line ({delimiters.size}, {n_delimiters_per_line}): {lines}", line_number=offending_line)
        # self._validated = True

    @classmethod
    def from_raw_buffer(cls, chunk: np.ndarray, header_data=None) -> "DelimitedBuffer":
        """Make EncodedArray of the chunk and extract all complete lines

        Also find all delimiters in the buffer

        Parameters
        ----------
        chunk : np.ndarray
            Raw bytes chunk as array
        header_data : 6
            Any header data that was read in `read_header`

        Returns
        -------
        DelimitedBuffer
            DelimitedBuffer object of all complete lines

        """
        chunk = EncodedArray(chunk, BaseEncoding)
        new_lines = np.flatnonzero(chunk == '\n')
        return cls(chunk[:new_lines[-1]+1], new_lines[:-1], header_data)
