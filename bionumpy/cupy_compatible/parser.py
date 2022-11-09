import cupy as cp

from ..io.parser import NumpyFileReader

import logging
logger = logging.getLogger(__name__)


class CupyFileReader(NumpyFileReader):
    def _get_buffer(self):
        b, size = super()._get_buffer()
        return cp.asanyarray(b), size
