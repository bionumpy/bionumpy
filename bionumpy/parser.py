import numpy as np
from .npdataclassstream import NpDataclassStream
import logging
from npstructures import npdataclass
from .sequences import EncodedArray
logger = logging.getLogger(__name__)


def wrapper(x):
    return x


def repr_bytes(n):
    if n < 10 ** 4:
        return str(n) + "b"
    elif n < 10 ** 7:
        return str(n // 1000) + "kb"
    elif n < 10 ** 11:
        return str(n // 1000000) + "Mb"
    return str(n // 1000000000) + "Gb"


class NumpyFileReader:
    """
    Class that handles the main task of reading in chunks of data
    from file that corresponds to full entries.
    """
    def __init__(self, file_obj, buffer_type, has_header=False):

        self._file_obj = file_obj
        """Inits the file reader with an opened file object, and a class of the buffer type used to parse it

        Parameters
        ----------
        file_obj : file
            Opened file object
        buffer_type : FileBuffer
            Buffer class used to parse the data
        has_header : bool
            whether or not the file has a header

        Examples
        --------
        FIXME: Add docs.

        """
        self._is_finished = False
        self._buffer_type = buffer_type
        self._has_header = has_header
        self._f_name = (
            self._file_obj.name
            if hasattr(self._file_obj, "name")
            else str(self._file_obj)
        )
        self._remove_initial_comments()
        self._header_data = self._buffer_type.read_header(self._file_obj)
        self._total_bytes_read = 0
        self._do_prepend = False
        self._prepend = []

    def __enter__(self):
        return self

    def __exit__(self):
        self._file_obj.close()

    def __iter__(self):
        return self.read_chunks()

    def set_prepend_mode(self):
        self._do_prepend = True

    def read(self):
        # self._remove_initial_comments()
        # self._header_data = self._buffer_type.read_header(self._file_obj)
        chunk = np.frombuffer(self._file_obj.read(), dtype=np.uint8)
        chunk, _  = self.__add_newline_to_end(chunk, chunk.size)
        return self._buffer_type.from_raw_buffer(chunk, header_data=self._header_data)

    def read_chunk(self, chunk_size=5000000):
        chunk = self.__get_buffer(chunk_size)
        if chunk is None:
            return None
        if len(self._prepend):
            chunk = np.concatenate((self._prepend, chunk))
            self._prepend = []

        self._total_bytes_read += chunk.size
        buff = self._buffer_type.from_raw_buffer(chunk, header_data=self._header_data)
        if not self._is_finished:
            if not self._do_prepend:
                self._file_obj.seek(buff.size - chunk.size, 1)
            else:
                self._prepend = chunk[buff.size:]
        if chunk is not None and chunk.size:
            return wrapper(buff)

    def read_chunks(self, chunk_size=5000000):
        #self._remove_initial_comments()
        #self._header_data = self._buffer_type.read_header(self._file_obj)
        while not self._is_finished:
            yield self.read_chunk(chunk_size)

    def __add_newline_to_end(self, chunk, bytes_read):
        if chunk[bytes_read - 1] != "\n":
            chunk = np.append(chunk, np.uint8(ord("\n")))
            bytes_read += 1
        if hasattr(self._buffer_type, "_new_entry_marker"):
            chunk = np.append(chunk, np.uint8(ord(self._buffer_type._new_entry_marker)))
            bytes_read += 1
        return chunk, bytes_read

    def __get_buffer(self, chunk_size):
        a, bytes_read = self.__read_raw_chunk(chunk_size)
        self._is_finished = bytes_read < chunk_size
        if bytes_read == 0:
            return None

        # Ensure that the last entry ends with newline. Makes logic easier later
        if self._is_finished:
            a, bytes_read = self.__add_newline_to_end(a, bytes_read)
        # 
        #     if a[bytes_read - 1] != ord("\n"):
        #         a = np.append(a, ord("\n"))
        #         bytes_read += 1
        #     if hasattr(self._buffer_type, "_new_entry_marker"):
        #         a = np.append(a, self._buffer_type._new_entry_marker)
        #         bytes_read += 1
        return a[:bytes_read]

    def __read_raw_chunk(self, chunk_size):
        b = np.frombuffer(self._file_obj.read(chunk_size), dtype="uint8")
        # assert not np.any(b & np.uint8(128)), "Unicdoe byte detected, not currently supported"
        return b.view(EncodedArray), b.size

    def _remove_initial_comments(self):
        if self._buffer_type.COMMENT == 0:
            return
        comment = self._buffer_type.COMMENT
        if isinstance(comment, str):
            comment = ord(comment)
        for line in self._file_obj:
            if line[0] != comment:
                self._file_obj.seek(-len(line), 1)
                break


class NpBufferedWriter:
    """ 
    File writer that can write @npdataclass objects
    to file
    """

    def __init__(self, file_obj, buffer_type):
        self._file_obj = file_obj
        self._buffer_type = buffer_type
        self._f_name = (
            self._file_obj.name
            if hasattr(self._file_obj, "name")
            else str(self._file_obj)
        )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._file_obj:
            self._file_obj.close()

    def close(self):
        self._file_obj.close()

    def write(self, data: npdataclass):
        """Write the provided data to file

        Parameters
        ----------
        data : npdataclass
            Data set containing entries

        """
        if isinstance(data, NpDataclassStream):
            for buf in data:
                self.write(buf)
            return 
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
