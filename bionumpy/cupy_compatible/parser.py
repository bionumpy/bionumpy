import sys
import cupy as cp

from ..io.parser import NumpyFileReader

import logging
logger = logging.getLogger(__name__)


class CupyFileReader(NumpyFileReader):
    def _get_buffer(self, min_chunk_size=5000000, max_chunk_size=None):
        temp = super()._get_buffer(min_chunk_size, max_chunk_size)
        if temp is None:
            return temp

        b = cp.asanyarray(temp)
        return b
