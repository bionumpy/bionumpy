import numpy as np
import logging
logger = logging.getLogger(__name__)

def repr_bytes(n):
    if n< 10**4:
        return str(n)
    elif n<10**7:
        return str(n//1000) + "kb"
    elif n<10**11:
        return str(n//1000000) + "Mb"
    return str(n//1000000000) + "Gb"

class NpBufferStream:
    def __init__(self, file_obj, buffer_type, chunk_size=5000000):
        self._file_obj = file_obj
        self._chunk_size = chunk_size
        self._is_finished = False
        self._buffer_type = buffer_type
        self._f_name = self._file_obj.name if hasattr(self._file_obj, "name") else str(self._file_obj)

    def __get_buffer(self):
        a, bytes_read = self.__read_raw_chunk()
        self._is_finished = bytes_read < self._chunk_size
        if bytes_read == 0:
            return None
        
        # Ensure that the last entry ends with newline. Makes logic easier later
        if self._is_finished  and a[bytes_read-1] != ord("\n"):
            a[bytes_read] = ord("\n")
            bytes_read += 1
        return a[:bytes_read]

    def __read_raw_chunk(self):
        array = np.empty(self._chunk_size, dtype="uint8")
        bytes_read = self._file_obj.readinto(array)
        return array, bytes_read

    def remove_initial_comments(self):
        if self._buffer_type.COMMENT == 0:
            return 
        for line in self._file_obj:
            if line[0] != self._buffer_type.COMMENT:
                self._file_obj.seek(-len(line), 1)
                break

    def __iter__(self):
        self.remove_initial_comments()
        chunk = self.__get_buffer()
        total_bytes = 0

        while not self._is_finished:
            total_bytes += chunk.size
            logger.debug(f"Read chunk of size {repr_bytes(chunk.size)} from {self._f_name}. (Total {repr_bytes(total_bytes)})")
            buff = self._buffer_type.from_raw_buffer(chunk)
            self._file_obj.seek(buff.size-self._chunk_size, 1)
            yield buff
            chunk = self.__get_buffer()
        if chunk is not None and chunk.size:
            yield self._buffer_type.from_raw_buffer(chunk)

class NpBufferedWriter:
    def __init__(self, file_obj, buffer_type):
        self._file_obj = file_obj
        self._buffer_type = buffer_type
        self._f_name = self._file_obj.name if hasattr(self._file_obj, "name") else str(self._file_obj)

    def write(self, data):
        # TODO: this is ugly
        bytes_array = self._buffer_type.from_data(data)
        self._file_obj.write(bytes(bytes_array))# .tofile(self._file_obj)
        logger.debug(f"Wrote chunk of size {repr_bytes(bytes_array.size)} to {self._f_name}")
