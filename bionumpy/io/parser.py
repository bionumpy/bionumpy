import logging
import numpy as np
from typing.io import IO
from npstructures import npdataclass
from ..streams import BnpStream
from ..encoded_array import EncodedArray
from ..streams.grouped import grouped_stream
from .file_buffers import FileBuffer
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
    def __init__(self, file_obj: IO, buffer_type: FileBuffer, has_header: bool = False):
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

        self._file_obj = file_obj
        self._is_finished = False
        self._buffer_type = buffer_type
        self._has_header = has_header
        self._f_name = (
            self._file_obj.name
            if hasattr(self._file_obj, "name")
            else str(self._file_obj)
        )
        self._header_data = self._buffer_type.read_header(self._file_obj)
        self._total_bytes_read = 0
        self._do_prepend = False
        self._prepend = []
        self.n_bytes_read = 0
        self.n_lines_read = 0


    def __enter__(self):
        return self

    def __exit__(self):
        self._file_obj.close()

    def __iter__(self):
        return self.read_chunks()

    def set_prepend_mode(self):
        self._do_prepend = True

    def read(self):
        chunk = np.frombuffer(self._file_obj.read(), dtype=np.uint8)
        chunk, _  = self.__add_newline_to_end(chunk, chunk.size)
        return self._buffer_type.from_raw_buffer(chunk, header_data=self._header_data)

    def read_chunk(self, min_chunk_size: int = 5000000, max_chunk_size: int = None) -> FileBuffer:
        """Read buffers of size `min_chunk_size` until at least one complete entry is found` Stop at `max_chunk_size`

        Read buffers from the file object. If the `from_raw_buffer`
        method of `self._file_obj` finds a complete entry, then read
        some more

        Parameters
        ----------
        min_chunk_size : int
            The size of the buffers to read
        max_chunk_size : int
            The max size of the combined buffers. If this is reached
            before a complete entry is found, raise Exception

        Returns
        -------
        np.ndarray
            A FileBuffer, containing at least one complete entry

        Examples
        --------
        FIXME: Add docs.

        """
        complete_entry_found = False
        temp_chunks = []
        local_bytes_read = 0
        local_lines_read = 0
        if len(self._prepend):
            temp_chunks.append(self._prepend)
        while not complete_entry_found:
            chunk = self._get_buffer(min_chunk_size, max_chunk_size)
            if chunk is None:
                return None
            temp_chunks.append(chunk)
            if max_chunk_size is not None and sum(chunk.size for chunk in chunks) > max_chunk_size:
                raise Exception("No complete entry found")
            local_bytes_read += chunk.size
            #local_lines_read += chunk.n_lines
            complete_entry_found = self._buffer_type.contains_complete_entry(temp_chunks)
            
        chunk = np.concatenate(temp_chunks)
        buff = self._buffer_type.from_raw_buffer(chunk, header_data=self._header_data)
        self._prepend = []
        if not self._is_finished:
            if not self._do_prepend:
                self._file_obj.seek(buff.size - chunk.size, 1)
            else:
                self._prepend = chunk[buff.size:]
        if chunk is not None and chunk.size:
            self.n_bytes_read += buff.size
            self.n_lines_read += buff.n_lines
            return wrapper(buff)

    def read_chunks(self, min_chunk_size: int = 5000000, max_chunk_size: int = None):
        while not self._is_finished:
            chunk = self.read_chunk(min_chunk_size, max_chunk_size)
            if chunk is None:
                break
            yield chunk

    def close(self):
        self._file_obj.close()

    def __add_newline_to_end(self, chunk, bytes_read):
        if chunk[bytes_read - 1] != ord("\n"):
            chunk = np.append(chunk, np.uint8(ord("\n")))
            bytes_read += 1
        if hasattr(self._buffer_type, "_new_entry_marker"):
            chunk = np.append(chunk, np.uint8(ord(self._buffer_type._new_entry_marker)))
            bytes_read += 1
        return chunk, bytes_read

    def _get_buffer(self, min_chunk_size: int = 5000000, max_chunk_size: int = None):
        a, bytes_read = self.__read_raw_chunk(min_chunk_size, max_chunk_size)
        self._is_finished = bytes_read < min_chunk_size
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

    def __read_raw_chunk(self, min_chunk_size: int = 5000000, max_chunk_size: int = None):
        b = np.frombuffer(self._file_obj.read(min_chunk_size), dtype="uint8")
        # assert not np.any(b & np.uint8(128)), "Unicdoe byte detected, not currently supported"
        return b, b.size


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
            dataset containing entries

        """
        if hasattr(self._buffer_type, 'make_header') and (not hasattr(self._file_obj, "mode") or self._file_obj.mode != 'ab') and getattr(self._buffer_type, 'INCLUDE_HEADER', False):
            header_array = self._buffer_type.make_header(data)
            self._file_obj.write(header_array)

        if isinstance(data, BnpStream):
            for buf in data:
                if len(buf) > 0:
                    self.write(buf)
            return
        if isinstance(data, grouped_stream):
            for name, buf in data:
                if len(buf) > 0:
                    self.write(buf)
            return
        bytes_array = self._buffer_type.from_data(data)
        if isinstance(bytes_array, EncodedArray):
            bytes_array = bytes_array.raw()
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
