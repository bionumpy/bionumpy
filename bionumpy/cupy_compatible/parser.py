import numpy as np
import cupy as cp

from ..parser import NumpyFileReader

import logging

logger = logging.getLogger(__name__)
wrapper = lambda x: x

def repr_bytes(n):
    if n < 10 ** 4:
        return str(n) + "b"
    elif n < 10 ** 7:
        return str(n // 1000) + "kb"
    elif n < 10 ** 11:
        return str(n // 1000000) + "Mb"
    return str(n // 1000000000) + "Gb"

class CpBufferStream(NumpyFileReader):

    def __iter__(self):
        self._remove_initial_comments()
        self._remove_header()
        chunk = self.__get_buffer()
        total_bytes = 0

        while not self._is_finished:
            total_bytes += chunk.size
            logger.debug(
                f"Read chunk of size {repr_bytes(chunk.size)} from {self._f_name}. (Total {repr_bytes(total_bytes)})"
            )
            buff = self._buffer_type.from_raw_buffer(chunk)
            self._file_obj.seek(buff.size - self._chunk_size, 1)
            yield wrapper(buff)
            chunk = self.__get_buffer()
        if chunk is not None and chunk.size:
            yield self._buffer_type.from_raw_buffer(chunk)

    def __get_buffer(self):
        a, bytes_read = self.__read_raw_chunk()
        self._is_finished = bytes_read < self._chunk_size
        if bytes_read == 0:
            return None

        # Ensure that the last entry ends with newline. Makes logic easier later
        if self._is_finished and a[bytes_read - 1] != ord("\n"):
            a = cp.append(a, ord("\n"))
            bytes_read += 1
        return a[:bytes_read]

    def __read_raw_chunk(self):
        b = np.frombuffer(self._file_obj.read(self._chunk_size), dtype="uint8")
        b = cp.asanyarray(b)
        return b, b.size
