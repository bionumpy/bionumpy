import numpy as np
my_empty =  np.empty

import logging

from npstructures import npdataclass

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


class NpBufferStream:
    """
    Class that handles the main task of reading in chunks of data
    from file that corresponds to full entries.
    """

    def __init__(self, file_obj, buffer_type, chunk_size=5000000):
        self._file_obj = file_obj
        self._chunk_size = chunk_size
        self._is_finished = False
        self._buffer_type = buffer_type
        self._f_name = (
            self._file_obj.name
            if hasattr(self._file_obj, "name")
            else str(self._file_obj)
        )

    def __iter__(self):
        self._remove_initial_comments()
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
            a[bytes_read] = ord("\n")
            bytes_read += 1
        return a[:bytes_read]

    def __read_raw_chunk(self):
        b = np.frombuffer(self._file_obj.read(self._chunk_size), dtype="uint8")
        return b, b.size
        # array = my_empty(self._chunk_size, dtype="uint8")
        # bytes_read = self._file_obj.readinto(array)
        #return array, bytes_read

    def _remove_initial_comments(self):
        if self._buffer_type.COMMENT == 0:
            return
        for line in self._file_obj:
            if line[0] != self._buffer_type.COMMENT:
                self._file_obj.seek(-len(line), 1)
                break


class NpBufferedWriter:
    def __init__(self, file_obj, buffer_type):
        self._file_obj = file_obj
        self._buffer_type = buffer_type
        self._f_name = (
            self._file_obj.name
            if hasattr(self._file_obj, "name")
            else str(self._file_obj)
        )

    def write(self, data: npdataclass):
        """Write the provided data to file

        Parameters
        ----------
        data : npdataclass
            Data set containing entries

        """
        bytes_array = self._buffer_type.from_data(data)
        self._file_obj.write(bytes(bytes_array))  # .tofile(self._file_obj)
        self._file_obj.flush()
        logger.debug(
            f"Wrote chunk of size {repr_bytes(bytes_array.size)} to {self._f_name}"
        )

def chunk_lines(stream, n_lines):
    cur_buffers = []
    remaining_lines = n_lines
    for chunk in stream:
        n_lines_in_chunk = len(chunk)
        while n_lines_in_chunk >= remaining_lines:
            cur_buffers.append(chunk[:remaining_lines])
            yield np.concatenate(cur_buffers)
            cur_buffers = []
            chunk = chunk[remaining_lines:]
            remaining_lines = n_lines
            n_lines_in_chunk = len(chunk)
        cur_buffers.append(chunk)
        remaining_lines -= n_lines_in_chunk
    yield np.concatenate(cur_buffers)
