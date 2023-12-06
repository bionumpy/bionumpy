from itertools import takewhile, repeat

from .delimited_buffers import DelimitedBufferWithInernalComments
from .parser import NumpyFileReader
from .. import config
from ..bnpdataclass import BNPDataClass
from .exceptions import FormatException
from ..datatypes import GTFEntry, VCFGenotypeEntry
from ..streams import NpDataclassStream
from ..bnpdataclass.lazybnpdataclass import create_lazy_class, ItemGetter


class NpDataclassReader:
    """
    This is the main file reader class for bionumpy. Ordinarily this should
    be created by using the `bnp.open` function that will create a reader with
    the correct attributes according to the file suffix. But this class can be used
    for instance if the a file name is not available (stdin), or you want more control over
    the reading.
    """
    def __init__(self, numpyfilereader: NumpyFileReader, lazy=None):
        self._reader = numpyfilereader
        self.__lazy_class = None
        self._lazy=lazy

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._reader.close()

    def close(self):
        self._reader.close()

    def read(self) -> BNPDataClass:
        """Read the whole file into a dataclass

        Use this for small files that can be held in memory

        Returns
        -------
        bnpdataclass
            A dataclass holdin all the entries in the class

        Examples
        --------
        4

        """
        chunk = self._reader.read()
        should_be_lazy = self._should_be_lazy(chunk)
        if should_be_lazy:
            return self._get_lazy_class(chunk.dataclass, chunk.header_data)(ItemGetter(chunk, chunk.dataclass))
        if chunk is None:
            return self._reader._buffer_type.dataclass.empty()
        return chunk.get_data()

    def _get_lazy_class(self, dataclass, header=None):
        if self.__lazy_class is None:
            self.__lazy_class = create_lazy_class(dataclass, header=header)
        return self.__lazy_class

    def _should_be_lazy(self, chunk):
        if ((not config.LAZY) and self._lazy is None) or (self._lazy is False):
            return False
        should_be_lazy = False
        if hasattr(chunk, 'get_field_by_number') and hasattr(chunk, 'dataclass'):
            if not issubclass(chunk.dataclass, (GTFEntry)):
                should_be_lazy = True
        return should_be_lazy

    def read_chunk(self, min_chunk_size: int = 5000000, max_chunk_size: int = None) -> BNPDataClass:
        """Read a single chunk into memory

        Read all complete entries in the next `chunk_size` bytes
        of the file. Useful for testing out algorithms on a small
        part of the file.

        Parameters
        ----------
        chunk_size: int
            How many bytes to read from file

        Returns
        -------
        bnpdataclass
            A dataclass holdin all the entries in the next chunk


        Examples
        --------
        5

        """
        n_lines_read = self._reader.n_lines_read
        chunk = self._reader.read_chunk(min_chunk_size, max_chunk_size)
        if chunk is None:
            return self._reader._buffer_type.dataclass.empty()
        try:
            if self._should_be_lazy(chunk):
                return self._get_lazy_class(chunk.dataclass, header=chunk.header_data)(ItemGetter(chunk, chunk.dataclass, n_lines_read))

            data = chunk.get_data()
        except FormatException as e:
            e.line_number += n_lines_read
            raise e
        return data

    def read_chunks(self, min_chunk_size: int = 5000000, max_chunk_size: int = None) -> NpDataclassStream:
        """Read the whole file in chunks

        This returns a generator yielding all the entries in the file 
        divided into chunks. Can be combined with functions decorated with
        `@streamable` to apply the function to each chunk in turn

        Parameters
        ----------
        chunk_size : int
            Number of bytes to read per chunk

        Returns
        -------
        NpDataclassStream
            4

        Examples
        --------
        5

        """
        data_stream = takewhile(len, (self.read_chunk(min_chunk_size, max_chunk_size) for _ in repeat(None)))

        return NpDataclassStream(data_stream, dataclass=self._reader._buffer_type.dataclass)


    def __iter__(self) -> NpDataclassStream:
        """Iteratate over chunks in the file

        Returns
        -------
        NpDataclassStream
            3

        """
        return self.read_chunks()
