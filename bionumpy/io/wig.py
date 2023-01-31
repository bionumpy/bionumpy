from .delimited_buffers import DelimitedBuffer
from .strops import split
from .exceptions import FormatException
import numpy as np


class WigBuffer(DelimitedBuffer):

    def _validate(self):
        chunk = self._data
        delimiters = self._delimiters[1:]
        n_delimiters_per_line = (
            next(i for i, d in enumerate(delimiters) if chunk[d] == '\n') + 1
        )
        self._n_cols = n_delimiters_per_line
        should_be_new_lines = chunk[delimiters[n_delimiters_per_line-1::n_delimiters_per_line]]
        if delimiters.size % n_delimiters_per_line != 0 or np.any(should_be_new_lines != "\n"):
            offending_line = np.flatnonzero(should_be_new_lines != "\n")[0]
            lines = split(self._data, '\n')
            raise FormatException(f"Irregular number of delimiters per line ({delimiters.size}, {n_delimiters_per_line}): {lines}", line_number=offending_line)
        self._validated = True
