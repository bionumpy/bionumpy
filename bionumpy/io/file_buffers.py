from functools import lru_cache
from typing import Optional, List, Union, Any, Type, Tuple

import numpy as np
from io import FileIO
from npstructures import RaggedView
from npstructures.raggedshape import RaggedView2

from .exceptions import FormatException
from .strops import str_to_int, str_to_int_with_missing, str_to_float, str_to_float_with_missing
from ..bnpdataclass import bnpdataclass, BNPDataClass
from ..encoded_array import EncodedArray, EncodedRaggedArray, Encoding, as_encoded_array
from ..encodings import BaseEncoding
from ..string_array import as_string_array
from ..typing import SequenceID
from ..util import is_subclass_or_instance

NEWLINE = "\n"


def move_intervals_to_digit_array(data, starts, ends, fill_value):
    if len(starts) == 0:
        return np.zeros_like(data, shape=((0, 0)))
    max_chars = np.max(ends - starts)
    view_starts = (ends - max_chars)
    indices = view_starts[..., None] + np.arange(max_chars)
    array = data[indices.ravel()]
    zeroed, _ = RaggedView(np.arange(starts.size) * max_chars,
                           max_chars - (ends - starts)).get_flat_indices()
    array[zeroed] = fill_value
    return array.reshape((-1, max_chars))


def move_intervals_to_right_padded_array(data, starts, ends, fill_value, stop_at=None):
    lens = ends - starts
    max_chars = np.max(lens)
    indices = np.minimum(starts[..., None] + np.arange(max_chars), data.size - 1)
    array = data[indices]
    del indices
    if stop_at is not None:
        new_lens = np.argmax(array == stop_at, axis=-1)
        lens = np.where(new_lens > 0, np.minimum(lens, new_lens), lens)
        max_chars = np.max(lens)
        array = array[:, :max_chars].ravel()
    z_lens = max_chars - lens
    row_idxs = np.flatnonzero(z_lens)
    if len(row_idxs) == 0:
        return array.reshape((-1, max_chars))
    first_row = row_idxs[0]
    cm = np.cumsum(z_lens[row_idxs])
    diffs = np.diff(row_idxs) * max_chars - z_lens[row_idxs[1:]]
    del row_idxs
    index_builder = np.ones(cm[-1], dtype=np.int32)
    index_builder[cm[:-1]] = diffs + 1
    del cm
    index_builder[0] = first_row * max_chars + lens[first_row]
    np.cumsum(index_builder, out=index_builder)
    zeroed = index_builder
    tmp = array.ravel()
    tmp[zeroed] = fill_value
    return array.reshape((-1, max_chars))


def wierd_padding(lens, max_chars):
    z_len = max_chars - lens
    mask = z_len != 0
    z_len = z_len[mask]
    cumsum = np.cumsum(z_len)
    L = cumsum[-1]
    zeroed = np.ones(L, dtype=int)
    row_idxs = np.flatnonzero(mask)
    diffs = np.diff(row_idxs) * max_chars
    diffs -= np.diff(z_len)
    zeroed[cumsum[:-1]] = diffs
    zeroed[0] = (row_idxs[0] + 1) * max_chars - z_len[0]
    zeroed = np.cumsum(zeroed)
    return zeroed


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
    supports_modified_write = True
    COMMENT = 0

    def __init__(self, data: EncodedArray, new_lines: np.ndarray):
        self._data = data
        self._new_lines = np.asanyarray(new_lines)
        self._is_validated = False

    @property
    @lru_cache()
    def size(self) -> int:
        return self.data.size

    @property
    def data(self) -> EncodedArray:
        return self._buffer_extractor.data

    def __getitem__(self, idx: Union[int, slice, List[int]]):
        return NotImplemented

    def raise_if(condition, *args, **kwargs):
        if condition:
            raise FormatException(*args, **kwargs)

    @property
    def header_data(self) -> Any:
        if hasattr(self, "_header_data"):
            return self._header_data
        return None

    @property
    def n_lines(self) -> int:
        return NotImplemented

    @classmethod
    def modify_class_with_header_data(cls, header_data: Any) -> Type["FileBuffer"]:
        return cls

    @classmethod
    def read_header(cls, file_object: FileIO) -> str:
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
        """Validate the buffer if it has not been validated yet"""
        if not self._is_validated:
            self._validate()
            self._is_validated = True

    def get_data(self) -> BNPDataClass:
        """Extract the data from the buffer

        The default way to extract data from the the buffer

        Returns
        -------
        BNPDataClass
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
        to_indices = ends[::-1, None] - max_chars + np.arange(max_chars)
        self._data[to_indices] = array[::-1]

    def _get_parser(self, field_type):
        parsers = [(str, lambda x: x),
                   (Encoding, lambda x: as_encoded_array(x, field_type)),
                   (SequenceID, as_string_array),
                   (int, lambda x: str_to_int(*x)),
                   (Optional[int], str_to_int_with_missing),
                   (bool, lambda x: str_to_int(x).astype(bool)),
                   (float, str_to_float),
                   (Optional[float], str_to_float_with_missing),
                   (List[int], lambda x: self._parse_split_ints(x)),
                   (List[float], lambda x: self._parse_split_floats(x)),
                   (List[bool], lambda x: self._parse_split_ints(x, sep="").astype(bool))]
        parser = None
        if is_subclass_or_instance(field_type, Encoding):
            parser = lambda x: as_encoded_array(x, field_type)
        for f, field_parser in parsers:
            if field_type == f:
                parser = field_parser
        return parser

    @classmethod
    def contains_complete_entry(cls, chunks: List[np.ndarray]) -> bool:
        n_new_lines = sum(np.count_nonzero(EncodedArray(chunk, BaseEncoding) == NEWLINE) for chunk in chunks)
        return n_new_lines >= cls.n_lines_per_entry

    @classmethod
    def process_field_for_write(cls, field_name, value):
        return value


class IncompleteEntryException(Exception):
    pass


class TextBufferExtractor:
    """
    Base class for extracting data from a text buffer.
    """
    def __init__(self, data: EncodedArray, field_starts: np.ndarray, field_ends: np.ndarray = None,
                 field_lens: np.ndarray = None):
        '''
        field_starts: n_entries x n_fields
        field_ends: n_entries x n_fields
        '''
        self._data = data
        self._field_starts = field_starts

        if field_lens is None:
            assert field_ends is not None
            self._field_lens = field_ends - field_starts
        else:
            assert field_ends is None
            self._field_lens = field_lens
        self._n_fields = field_starts.shape[1]

    @property
    def data(self) -> EncodedArray:
        return self._data

    @property
    def n_fields(self) -> int:
        return self._n_fields

    def __len__(self):
        return len(self._field_starts)

    def __getitem__(self, idx: Union[int, slice, List[int]]) -> 'TextBufferExtractor':
        return self.__class__(self._data,
                              field_starts=self._field_starts[idx],
                              field_lens=self._field_lens[idx])

    def get_field_by_number(self, field_nr: int, keep_sep: bool=False) -> EncodedRaggedArray:
        """
        Extract the data for a single field.
        Parameters
        ----------
        field_nr: int
        keep_sep: bool

        Returns
        -------
        EncodedRaggedArray

        """
        assert field_nr < self._n_fields, (field_nr, self._n_fields)
        lens = self._field_lens.ravel()[field_nr::self._n_fields]
        if keep_sep:
            lens = lens + 1
        starts = self._field_starts.ravel()[field_nr::self._n_fields]
        return self._extract_data(lens, starts)

    def _extract_data(self, lens, starts):
        values = EncodedRaggedArray(self._data, RaggedView2(starts, lens))
        assert len(values) == len(self), (self._field_starts, self._field_lens, self._n_fields, self._data)
        return values

    def get_fixed_length_field(self, field_nr: int, field_length: int)-> EncodedArray:
        indices = self._field_starts[:, field_nr, None] + np.arange(field_length)
        return self._data[indices]

    def get_padded_field(self, field_nr, stop_at=None) -> EncodedArray:
        starts = self._field_starts[:, field_nr]
        if starts.size == 0:
            return np.zeros_like(self._data, shape=(len(starts), 0))
        lens = self._field_lens[:, field_nr]
        ends = lens + starts

        array = move_intervals_to_right_padded_array(self._data, starts.ravel(), ends.ravel(), fill_value='\x00',
                                                     stop_at=stop_at)
        return array.reshape(starts.shape + (array.shape[-1],))

    def get_digit_array(self, field_nr: int) -> Tuple[EncodedArray, Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Extract the digits of the field as a 2D array of encoded integres.

        Parameters
        ----------
        field_nr: int

        Returns
        -------
        EncodedArray
        """

        starts = self._field_starts[:, field_nr]
        possible_signs = self._data[starts]
        is_negative = possible_signs == "-"
        is_positive = possible_signs == "+"
        if np.any(is_negative) or np.any(is_positive):
            return self.get_field_by_number(field_nr), is_negative, is_positive
        digit_array = move_intervals_to_digit_array(self._data, starts, starts + self._field_lens[:, field_nr],
                                                    fill_value='0')
        return digit_array, None, None

    @classmethod
    def concatenate(cls, buffers: List['TextBufferExtractor']):
        """
        Concatenate multiple buffers into a single buffer.

        Parameters
        ----------
        buffers: List[TextBufferExtractor]

        Returns
        -------
        TextBufferExtractor

        """
        sizes = np.array([b._data.size for b in buffers])
        offsets = np.insert(np.cumsum(sizes), 0, 0)
        data = np.concatenate([b._data for b in buffers])
        starts = np.concatenate([b._field_starts + offset for b, offset in zip(buffers, offsets)])
        lens = np.concatenate([b._field_lens for b in buffers])
        return cls(data, starts, field_lens=lens)


class TextThroughputExtractor(TextBufferExtractor):
    """
    TextBufferExtractor made especially for making it fast to write a modified or filtered buffer to file again.
    """

    def __init__(self, data: EncodedArray, field_starts: np.ndarray, field_ends: np.ndarray = None, field_lens=None,
                 entry_starts: np.ndarray = None, entry_ends: np.ndarray = None, is_contiguous=True):
        if field_lens is None:
            field_lens = field_ends - field_starts
        super().__init__(data, field_starts, field_lens=field_lens)
        self._entry_starts = entry_starts
        self._entry_ends = entry_ends
        self._is_contiguous = is_contiguous

    @classmethod
    def concatenate(cls, buffers: List['TextThroughputExtractor']):
        sizes = np.array([b._data.size for b in buffers])
        offsets = np.insert(np.cumsum(sizes), 0, 0)
        data = np.concatenate([b._data for b in buffers])
        starts = np.concatenate([b._field_starts + offset for b, offset in zip(buffers, offsets)])
        lens = np.concatenate([b._field_lens for b in buffers])
        entry_starts = np.concatenate([b._entry_starts + offset for b, offset in zip(buffers, offsets)])
        entry_ends = np.concatenate([b._entry_ends + offset for b, offset in zip(buffers, offsets)])
        return cls(data, starts, field_lens=lens, entry_starts=entry_starts, entry_ends=entry_ends,
                   is_contiguous=all(b._is_contiguous for b in buffers))

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
        self._field_starts = self._field_starts - offsets[:, None]
        self._is_contiguous = True

    @property
    def data(self) -> EncodedArray:
        if not self._is_contiguous:
            self._make_contigous()
        return self._data

    def get_fields_by_range(self, from_nr: Optional[int] = None, to_nr: Optional[int] = None, keep_sep=False):
        assert from_nr is not None
        assert to_nr is None
        starts = self._field_starts[:, from_nr]
        lens = self._entry_ends - starts
        if not keep_sep:
            lens -= 1
        return self._extract_data(lens, starts)
