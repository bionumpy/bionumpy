from functools import lru_cache

import numpy as np
from io import FileIO
from npstructures import RaggedView
from npstructures.raggedshape import RaggedView2

from .exceptions import FormatException
from ..bnpdataclass import bnpdataclass
from ..encoded_array import EncodedArray, EncodedRaggedArray
from ..encodings import BaseEncoding

NEWLINE = "\n"


def move_intervals_to_digit_array(data, starts, ends, fill_value):
    max_chars = np.max(ends - starts)
    view_starts = (ends - max_chars)
    indices = view_starts[:, None] + np.arange(max_chars)
    array = data[indices.ravel()]
    zeroed, _ = RaggedView(np.arange(starts.size) * max_chars, max_chars - (ends - starts)).get_flat_indices()
    array[zeroed] = fill_value
    return array.reshape((-1, max_chars))


class FileBuffer:
    """ Base class for file buffer classes. Different FileBuffer classes
    should correspond to different file formats.

    A FileBuffer class can extract the text/bytes corresponding to complete
    entries from a raw buffer (`from_raw_buffer`) and convert the bytes/text from
    those entries into data in `@npdataclass` objects.

    The `from_data` method should convert `@npdataclass` objects into text/bytes
    that can be written to file.

    This base class provides some convenience methods to extract text form parts of
    a buffer into meaningful data
    """

    _buffer_divisor = 1
    COMMENT = 0

    def __init__(self, data: EncodedArray, new_lines: np.ndarray):
        self._data = data
        self._new_lines = np.asanyarray(new_lines)
        self._is_validated = False

    @property
    @lru_cache()
    def size(self):
        return self.data.size

    @property
    def data(self):
        return self._buffer_extractor.data

    def __getitem__(self, idx):
        return NotImplemented

    def raise_if(condition, *args, **kwargs):
        if condition:
            raise FormatException(*args, **kwargs)

    @property
    def header_data(self):
        if hasattr(self, "_header_data"):
            return self._header_data
        return None

    @property
    def n_lines(self):
        return NotImplemented
        return len(self._new_lines)

    @classmethod
    def modify_class_with_header_data(cls, header_data):
        return cls
        # class NewClass(cls):
        #     _header_data = header_data


    @classmethod
    def read_header(cls, file_object: FileIO):
        """Read the header data from the file

        The data returned here is passed to each buffer through the
        `from_raw_buffer` method when reading a file, and can so influence
        how the data in the file is parsed.

        This function should leave the file pointer pointing to the beginning
        of the data to be read.

        Parameters
        ----------
        file_object : file
            The file object to read the file from

        Examples
        --------
        6

        """
        if cls.COMMENT == 0:
            return
        header = []
        comment = cls.COMMENT
        if isinstance(comment, str):
            comment = ord(comment)
        for line in file_object:
            if line[0] != comment:
                file_object.seek(-len(line), 1)
                break
            header.append(line.decode("utf-8"))
        return "".join(header)

    @classmethod
    def from_raw_buffer(cls, raw_buffer: np.ndarray, header_data=None) -> "FileBuffer":
        """Create a buffer with full entries

        A raw buffer can end with data that does not represent full entries.
        This method extracts all the full entries, so that the next buffer can
        start from the last incomplete entry.

        Parameters
        ----------
        chunk : np.ndarray
            Raw buffer with data that might end with incomplete entry

        Returns
        -------
        'FileBuffer'
            Buffer with complete entries

        Examples
        --------
        """

        return NotImplemented

    @classmethod
    def from_data(cls, data: bnpdataclass) -> "FileBuffer":
        """Create FileBuffer from a dataset

        Create a FileBuffer that can be written to file

        Parameters
        ----------
        data : npdataclass
            dataset containing the data to be written

        Returns
        -------
        'FileBuffer'
            FileBuffer containing the data
        """
        return NotImplemented

    def validate_if_not(self):
        if not self._is_validated:
            self._validate()
            self._is_validated = True

    def get_data(self) -> bnpdataclass:
        """Extract the data from the buffer

        The default way to extract data from the the buffer

        Returns
        -------
        npdataclass
            dataset containing the data from the buffer
        """
        return NotImplemented

    def _move_intervals_to_2d_array(self, starts, ends, fill_value=0):
        return move_intervals_to_digit_array(self._data, starts, ends, fill_value)

    def _move_intervals_to_ragged_array(self, starts, ends=None, lens=None, as_sequence=True):
        if lens is None:
            lens = ends - starts
        e = EncodedRaggedArray(self._data, RaggedView(starts, lens))
        e.is_contigous = False
        return e
        # return EncodedRaggedArray(self._data[indices], shape)

    def _move_2d_array_to_intervals(self, array, starts, ends):
        max_chars = array.shape[-1]
        to_indices = ends[::-1, None]-max_chars+np.arange(max_chars)
        self._data[to_indices] = array[::-1]

    @classmethod
    def contains_complete_entry(cls, chunks):
        n_new_lines = sum(np.count_nonzero(EncodedArray(chunk, BaseEncoding) == NEWLINE) for chunk in chunks)
        return n_new_lines >= cls.n_lines_per_entry


class IncompleteEntryException(Exception):
    pass


class TextBufferExtractor:
    def __init__(self, data: EncodedArray, field_starts: np.ndarray, field_ends: np.ndarray=None, field_lens: np.ndarray=None):
        '''
        field_starts: n_entries x n_fields
        field_ends: n_entries x n_fields
        '''
        self._data = data
        self._field_starts = field_starts

        if field_lens is None:
            assert field_ends is not None
            self._field_lens = field_ends-field_starts
        else:
            assert field_ends is None
            self._field_lens = field_lens
        self._n_fields = field_starts.shape[1]

    @property
    def data(self):
        return self._data

    def __len__(self):
        return len(self._field_starts)

    def __getitem__(self, idx):
        return self.__class__(self._data,
                              field_starts=self._field_starts[idx],
                              field_lens=self._field_lens[idx])

    def get_field_by_number(self, field_nr: int, keep_sep=False):
        assert field_nr < self._n_fields, (field_nr, self._n_fields)
        lens = self._field_lens.ravel()[field_nr::self._n_fields]
        if keep_sep:
            lens = lens + 1
        starts = self._field_starts.ravel()[field_nr::self._n_fields]
        values = EncodedRaggedArray(self._data, RaggedView2(starts, lens))
        # values = e[field_nr::self._n_fields]
        assert len(values) == len(self), (self._field_starts, self._field_lens, field_nr, self._n_fields, self._data)
        return values

    def get_fixed_length_field(self, field_nr: int, field_length: int):
        indices = self._field_starts[:, field_nr, None] + np.arange(field_length)
        return self._data[indices]

    def get_digit_array(self, field_nr: int):
        starts = self._field_starts[:, field_nr]
        possible_signs = self._data[starts]
        is_negative = possible_signs == "-"
        is_positive = possible_signs == "+"
        if np.any(is_negative) or np.any(is_positive):
            return self.get_field_by_number(field_nr), is_positive, is_negative
        digit_array = move_intervals_to_digit_array(self._data, starts, starts+self._field_lens[:, field_nr], fill_value='0')
        return digit_array, None, None


class TextThroughputExtractor(TextBufferExtractor):
    def __init__(self, data: EncodedArray, field_starts: np.ndarray, field_ends: np.ndarray=None, field_lens=None, entry_starts: np.ndarray=None, entry_ends:np.ndarray=None, is_contiguous=True):
        if field_lens is None:
            field_lens = field_ends-field_starts
        super().__init__(data, field_starts, field_lens=field_lens)
        self._entry_starts = entry_starts
        self._entry_ends = entry_ends
        self._is_contiguous = is_contiguous

    def __getitem__(self, idx):
        return self.__class__(self._data,
                              field_starts=self._field_starts[idx],
                              field_lens=self._field_lens[idx],
                              entry_starts=self._entry_starts[idx],
                              entry_ends=self._entry_ends[idx], is_contiguous=False)

    def _make_contigous(self):
        assert not self._is_contiguous
        lens = self._entry_ends - self._entry_starts
        new_starts = np.insert(np.cumsum(lens), 0, 0)
        offsets = self._entry_starts - new_starts[:-1]
        self._data = EncodedRaggedArray(self._data, RaggedView2(self._entry_starts, lens)).ravel()
        self._entry_starts = new_starts[:-1]
        self._entry_ends = new_starts[1:]
        self._field_starts = self._field_starts-offsets[:, None]

        # self._field_ends = self._field_ends-offsets[:, None]
        self._is_contiguous = True

    @property
    def data(self):
        if not self._is_contiguous:
            self._make_contigous()
        return self._data
        # return EncodedRaggedArray(self._data,
        #                           RaggedView2(self._entry_starts, self._entry_ends-self._entry_starts)).ravel()

