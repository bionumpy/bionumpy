import io
import logging
import dataclasses
from typing import List, Optional
from npstructures import RaggedArray, RaggedShape
from ..bnpdataclass import bnpdataclass, BNPDataClass
from ..bnpdataclass.lazybnpdataclass import LazyBNPDataClass
from ..datatypes import (Interval, SequenceEntry, Bed12, Bed6, BedGraph,
                         GTFEntry, GFFEntry, ChromosomeSize, NarrowPeak, GfaPath)
from ..encoded_array import EncodedArray, EncodedRaggedArray
from ..encoded_array import as_encoded_array
from ..encodings import Encoding
from ..encodings.exceptions import EncodingError
from ..encodings.alphabet_encoding import DigitEncoding
from ..encoded_array import BaseEncoding
from ..string_array import as_string_array
from ..typing import SequenceID
from ..util import is_subclass_or_instance
from .file_buffers import FileBuffer, NEWLINE, TextBufferExtractor, TextThroughputExtractor
from .strops import (
    split, str_to_int, str_to_float, str_to_int_with_missing, str_to_float_with_missing)
from .dump_csv import dump_csv, join_columns
from .exceptions import FormatException
import numpy as np
from ..bnpdataclass.bnpdataclass import make_dataclass
from ..bnpdataclass.lazybnpdataclass import create_lazy_class


class DelimitedBuffer(FileBuffer):
    """Base class for file buffers for delimited files such as csv or tsv.
    Each line should correspond to an entry, and each column to a variable.

    Provides convenience methods for extracting and decoding integers from columns,
    and text from columns into Sequences objects
    """

    DELIMITER = "\t"
    COMMENT = "#"
    HAS_UNCOMMENTED_HEADER_LINE = False
    n_lines_per_entry = 1

    def __init__(self, buffer_extractor: TextBufferExtractor, header_data=None):
        self._buffer_extractor = buffer_extractor
        self._header_data = header_data
        self._is_validated = True

    def concatenate(self, buffers):
        return self.__class__(self._buffer_extractor.concatenate([b._buffer_extractor for b in buffers]),
                              header_data=self._header_data)

    @classmethod
    def from_raw_buffer(cls, chunk: np.ndarray, header_data=None) -> "DelimitedBuffer":
        """Make EncodedArray of the chunk and extract all complete lines

        Also find all delimiters in the buffer

        Parameters
        ----------
        chunk : np.ndarray
            Raw bytes chunk as array
        header_data : 6
            Any header data that was read in `read_header`

        Returns
        -------
        DelimitedBuffer
            DelimitedBuffer object of all complete lines

        """
        chunk = EncodedArray(chunk, BaseEncoding)
        mask = chunk == NEWLINE
        mask |= chunk == cls.DELIMITER
        delimiters = np.flatnonzero(mask)
        entry_ends = np.flatnonzero(chunk[delimiters] == '\n')
        if entry_ends.size == 0:
            logging.warning("Foud no new lines. Chunk size may be too low. Try increasing")
            raise
        n_fields = cls._get_n_fields(entry_ends)
        size = delimiters[entry_ends[-1]] + 1
        delimiters = np.insert(delimiters[:entry_ends[-1] + 1], 0, -1)
        buffer_extractor = cls._get_buffer_extractor(
            chunk[:size], delimiters, n_fields)
        return cls(buffer_extractor, header_data)
        # return cls(chunk[:new_lines[-1] + 1], new_lines, delimiters, header_data, buffer_extractor=buffer_extractor)

    def __getitem__(self, idx):
        return self.__class__(self._buffer_extractor[idx], self._header_data)

    @property
    def entries(self):
        if not hasattr(self, "_entries"):
            lengths = np.diff(self._new_lines)
            lengths = np.insert(lengths, 0, self._new_lines[0] + 1)
            self._entries = EncodedRaggedArray(self._data, RaggedShape(lengths))
        return self._entries

    def get_text(self, col, fixed_length=True, keep_sep=False):
        assert not fixed_length and not keep_sep
        if not fixed_length and not keep_sep:
            return self._buffer_extractor.get_field_by_number(col)
        """Extract text from a column

        Extract strings from the specified column into either a 2d
        array or a RaggedArray

        Parameters
        ----------
        col : int
            column index
        fixed_length : bool
            whether all strings have equal length
        keep_sep : bool
            keep seperator at end

        Examples
        --------
        FIXME: Add docs.

        """

    @classmethod
    def join_fields(cls, fields_list: List[EncodedRaggedArray]) -> EncodedRaggedArray:
        return join_columns(fields_list, cls.DELIMITER).ravel()

    def get_field_range_as_text(self, *args, **kwargs) -> EncodedRaggedArray:
        return self.get_column_range_as_text(*args, **kwargs)

    def get_column_range_as_text(self, col_start, col_end, keep_sep=False) -> EncodedRaggedArray:
        """Get multiple columns as text

       Parameters
        ----------
        col_start : int
            column index to start at
        col_end : int
            column index to end at
        keep_sep : bool
            keep seperator at end

        Returns
        -------
        EncodedRaggedArray
            EncodedRaggedArray of the columns
        """
        self.validate_if_not()
        assert col_end == col_start + 1
        return self._buffer_extractor.get_field_by_number(col_start)

    @classmethod
    def from_data(cls, data: BNPDataClass) -> "DelimitedBuffer":
        """Put each field of the dataclass into a column in a buffer.

        Parameters
        data : bnpdataclass
            Data
        """
        if isinstance(data, LazyBNPDataClass):
            return cls.from_data(data.get_data_object())
        data_dict = [(field.type, getattr(data, field.name)) for field in dataclasses.fields(data)]
        return dump_csv(data_dict, cls.DELIMITER)

    @classmethod
    def make_header(cls, data: BNPDataClass):
        header = ""
        if data.has_context("header"):
            header = data.get_context("header")
        return bytes(header, "ascii")

    def get_data(self) -> BNPDataClass:
        """Parse the data in the buffer according to the fields in _dataclass

        Returns
        -------
        _dataclass
            Dataclass with parsed data

        """
        self.validate_if_not()
        columns = {}
        fields = dataclasses.fields(self.actual_dataclass)
        for col_number, field in enumerate(fields):
            col = self._get_field_by_number(col_number, field.type)
            columns[field.name] = col
        n_entries = len(next(col for col in columns if col is not None))
        columns = {c: value if c is not None else np.empty((n_entries, 0))
                   for c, value in columns.items()}
        data = self.actual_dataclass(**columns)
        data.set_context("header", self._header_data)
        return data

    @property
    def actual_dataclass(self):
        return self.dataclass

    def get_field_by_number(self, field_nr: int, field_type: type = object):
        """Get a field by number"""
        self.validate_if_not()
        if field_type is None:
            field_type = dataclasses.fields(self.actual_dataclass)[field_nr]
        return self._get_field_by_number(
            field_nr, field_type)

    def _get_field_by_number(self, col_number, field_type):

        if field_type is None:
            return None
        self.validate_if_not()
        if field_type == int:
            subresult = self._buffer_extractor.get_digit_array(col_number)
            text = subresult[0]
        elif field_type == SequenceID:
            subresult = self._buffer_extractor.get_padded_field(col_number)
            text = subresult
        else:
            subresult: EncodedRaggedArray = self._buffer_extractor.get_field_by_number(
                col_number,
                keep_sep=(field_type == List[int] or field_type == List[float]))
            text = subresult
        assert isinstance(text, (EncodedRaggedArray, EncodedArray)), text
        parser = self._get_parser(field_type)
        if parser is None:
            assert False, (self.__class__, field_type)
        try:
            parsed = parser(subresult)
            assert len(parsed) == len(text)
        except EncodingError as e:
            if isinstance(text, EncodedArray):
                row_number = e.offset // text.shape[1]
            else:
                row_number = np.searchsorted(np.cumsum(text.lengths), e.offset, side="right")
            raise FormatException(e.args[0], line_number=row_number)
        return parsed

    def count_entries(self) -> int:
        """Count the number of entries in the buffer"""
        return len(self._buffer_extractor)

    @property
    def n_lines(self) -> int:
        return len(self._buffer_extractor)

    def _parse_split_floats(self, text, sep=','):
        function = str_to_float
        return self._parse_split_fields(text, function, sep)

    def _parse_split_ints(self, text, sep=','):
        function = str_to_int
        return self._parse_split_fields(text, function, sep)

    def _parse_split_fields(self, text, function, sep):
        if len(sep):
            try:
                text[:, -1] = sep
            except ValueError:
                text = text.copy()
                text[:, -1] = sep
            int_strings = split(text.ravel()[:-1], sep=sep)

            if np.any(int_strings.lengths == 0):
                mask = int_strings.lengths != 0
                return RaggedArray(function(int_strings[mask]), (text == sep).sum(axis=-1),
                                   safe_mode=False)  # TODO: is it necessary with unsafe mode here
            return RaggedArray(function(int_strings), (text == sep).sum(axis=-1))
        else:
            mask = as_encoded_array(text.ravel(), DigitEncoding).raw()
            return RaggedArray(mask, text.shape)

    @classmethod
    def _get_n_fields(cls, entry_ends):
        return entry_ends[0] + 1

    @property
    def __buffer_extractor(self):
        if self.__buffer_extractor is None:
            self.__buffer_extractor = self._get_buffer_extractor()
        return self.__buffer_extractor

    @classmethod
    def _get_buffer_extractor(cls, data, delimiters, n_cols) -> TextThroughputExtractor:
        starts = delimiters[:-1].reshape(-1, n_cols) + 1
        ends = delimiters[1:].reshape(-1, n_cols)
        ends = cls._modify_for_carriage_return(ends, data)
        entry_starts = starts[:, 0]
        entry_ends = ends[:, -1] + 1
        return TextThroughputExtractor(data, starts, field_ends=ends, entry_starts=entry_starts, entry_ends=entry_ends)

    @classmethod
    def _modify_for_carriage_return(cls, ends, data):
        if data.size == 0 or ends[0, -1] == 0:
            return ends
        if data[ends[0, -1] - 1] == '\r':
            ends = ends.copy()
            ends[:, -1] -= data[ends[:, -1] - 1] == '\r'
        return ends

    @staticmethod
    def _move_ints_to_digit_array(ints, n_digits):
        powers = np.uint8(10) ** np.arange(n_digits)[::-1]
        ret = (ints[..., None] // powers) % 10
        return EncodedArray(ret, DigitEncoding)

    def _validate(self):
        chunk = self._data
        delimiters = self._delimiters[1:]
        n_delimiters_per_line = (
                next(i for i, d in enumerate(delimiters) if chunk[d] == NEWLINE) + 1
        )
        self._n_cols = n_delimiters_per_line
        should_be_new_lines = chunk[delimiters[n_delimiters_per_line - 1::n_delimiters_per_line]]
        if delimiters.size % n_delimiters_per_line != 0 or np.any(should_be_new_lines != "\n"):
            offending_line = np.flatnonzero(should_be_new_lines != "\n")[0]
            lines = split(self._data, '\n')
            raise FormatException(
                f"Irregular number of delimiters per line ({delimiters.size}, {n_delimiters_per_line}): {lines}",
                line_number=offending_line)
        self._validated = True


class GfaSequenceBuffer(DelimitedBuffer):
    dataclass = SequenceEntry
    def get_data(self):
        ids = self.get_text(1, fixed_length=False)
        sequences = self.get_text(col=2, fixed_length=False)
        return SequenceEntry(ids, sequences)

    def get_field_by_number(self, field_nr: int, field_type: type = object):
        return super().get_field_by_number(field_nr + 1, field_type)

    @classmethod
    def from_data(cls, data: SequenceEntry) -> EncodedArray:
        return dump_csv([(str, as_encoded_array(["S"] * len(data))),
                         (str, data.name),
                         (str, data.sequence)])


class GfaPathBuffer(DelimitedBuffer):
    def get_data(self):
        name = self.get_text(1, fixed_length=False)
        nodes_lists = self.get_text(2, keep_sep=True, fixed_length=False)
        nodes_lists[:, -1] = ","
        lengths = np.sum(nodes_lists == ",", axis=-1)
        all_node_texts = split(nodes_lists.ravel()[:-1], ",")
        int_text = all_node_texts[:, :-1]
        node_ids = str_to_int(int_text)
        directions = np.where(all_node_texts[:, -1] == "+", 1, -1)
        node_ids = RaggedArray(node_ids, lengths)
        directions = RaggedArray(directions, lengths)
        data = GfaPath(name, node_ids, directions)
        return data


def get_bufferclass_for_datatype(_dataclass: bnpdataclass, delimiter: str = "\t", has_header: bool = False,
                                 comment: str = "#",
                                 sub_delimiter=",") -> type:
    """Create a FileBuffer class that can read a delimited file with the fields specified in `_dataclass`

    This can be used to create a parser for a custom delimited file format and also more generic csv 
    reading. The order of the fields in the `_dataclass` is used as the order of the columns in the delimited
    file, unless `has_header=True`, in whcih case the name of the field corresponds to the name in the header.

    Parameters
    ----------
    _dataclass : bnpdataclass
        The dataclass used as template for the DelimitedBuffer
    delimiter : str
        The character used to separate the columns
    has_header : bool
        Wheter a header line should used to match the dataclass fields to columns
    comment : str
        The characted used to specify comment/unused lines in the file

    """

    class DatatypeBuffer(DelimitedBuffer):
        DELIMITER = delimiter
        COMMENT = comment
        HAS_UNCOMMENTED_HEADER_LINE = has_header
        dataclass = _dataclass

        # fields = None
        # data: EncodedArray, new_lines: np.ndarray, delimiters: np.ndarray = None,
        # def __init__(self, buffer_extractor, header_data: List[str] = None):
        #     super().__init__(buffer_extractor, header_data)
        #     # super().__init__(data, new_lines, delimiters, header_data)
        #   # self.set_fields_from_header(header_data)

        @classmethod
        def modify_class_with_header_data(cls, columns):
            fields = dataclasses.fields(cls.dataclass)
            type_dict = {field.name: field.type for field in fields}
            new_fields = [(name, type_dict[name]) if name in type_dict else (name, str) for name in columns]
            # ordered_fields = [next(field for field in fields if field.name == col) for col in columns]
            tmp = make_dataclass(
                new_fields, cls.dataclass.__name__ + 'Permuted')
            new_dataclass = make_dataclass([], cls.dataclass.__name__ + 'Permuted', bases=(cls.dataclass, tmp))
            assert [f.name for f in dataclasses.fields(tmp)] == columns, (
                columns, [f.name for f in dataclasses.fields(tmp)])

            class NewClass(cls):
                _actual_dataclass = cls.dataclass
                dataclass = tmp
                lazy_class = create_lazy_class(tmp)

            return NewClass

        def get_data(self) -> BNPDataClass:
            return super().get_data().astype(self._actual_dataclass)

        @classmethod
        def read_header(cls, file_object: io.FileIO) -> List[str]:
            """Read the column names from the header if `has_header=True`

            Parameters
            ----------
            file_object : io.FileIO

            Returns
            -------
            List[str]
                Column names
            """
            super().read_header(file_object)
            if not has_header:
                return None
            delimiter = cls.DELIMITER
            if not isinstance(delimiter, str):
                delimiter = chr(delimiter)
            return file_object.readline().decode('ascii').strip().split(delimiter)

        @classmethod
        def make_header(cls, data: bnpdataclass):
            """makes a header from field names separated by delimiter"""
            return bytes(cls.DELIMITER.join([field.name for field in dataclasses.fields(data)]) + "\n", 'ascii')

    DatatypeBuffer.__name__ = _dataclass.__name__ + "Buffer"
    DatatypeBuffer.__qualname__ = _dataclass.__qualname__ + "Buffer"
    return DatatypeBuffer


class BedBuffer(DelimitedBuffer):
    dataclass = Interval

    def __get_integers(self, cols: list) -> np.ndarray:
        ''' This is maybe a quicker way to parse ints than the default'''
        """Get integers from integer string

        Extract integers from the specified columns

        Parameters
        ----------
        cols : list
            list of columns containing integers

        """
        assert np.all(cols < self._n_cols), (str(self._data), cols, self._n_cols)
        cols = np.asanyarray(cols)
        assert cols.size == 1
        integer_starts = self._col_starts(cols)
        integer_ends = self._col_ends(cols)
        array = self._move_intervals_to_2d_array(integer_starts, integer_ends, fill_value='0')
        try:
            digits = as_encoded_array(array, DigitEncoding).raw()
        except EncodingError as e:
            row_number = e.offset // array.shape[-1]  # rows._shape.starts, e.offset, side="right")-1
            raise FormatException(e.args[0], line_number=row_number)
        powers = 10 ** np.arange(digits.shape[-1])[::-1]
        return digits.dot(powers).reshape(-1, cols.size)


class Bed6Buffer(BedBuffer):
    dataclass = Bed6


class Bed12Buffer(Bed6Buffer):
    dataclass = Bed12


class BdgBuffer(BedBuffer):
    dataclass = BedGraph


class NarrowPeakBuffer(DelimitedBuffer):
    dataclass = NarrowPeak


class GTFBuffer(DelimitedBuffer):
    dataclass = GTFEntry


class ChromosomeSizeBuffer(DelimitedBuffer):
    dataclass = ChromosomeSize


class DelimitedBufferWithInernalComments(DelimitedBuffer):

    @classmethod
    def _calculate_col_starts_and_ends(cls, data, delimiters):
        comment_mask = (data[delimiters[:-1]] == '\n') & (data[delimiters[:-1] + 1] == cls.COMMENT)
        comment_mask = np.flatnonzero(comment_mask)
        start_delimiters = np.delete(delimiters, comment_mask)[:-1]
        end_delimiters = np.delete(delimiters, comment_mask + 1)
        if data[0] != cls.COMMENT:
            start_delimiters = np.insert(start_delimiters, 0, -1)
        else:
            end_delimiters = end_delimiters[1:]
        return start_delimiters + 1, end_delimiters

    def _col_ends(self, col):
        return self._wig_col_ends[:, col]

    def _col_starts(self, col):
        return self._wig_col_starts[:, col]

    def _validate(self):
        self._validated = True

    def __init__(self, buffer_extractor: TextBufferExtractor, header_data=None):
        self._buffer_extractor = buffer_extractor
        self._header_data = header_data
        self._is_validated = True

        # data: EncodedArray, new_lines: np.ndarray, delimiters: np.ndarray = None, header_data=None):
        # delimiters_mask = data == self.DELIMITER
        # delimiters_mask[new_lines] = True
        # delimiters = np.append(np.flatnonzero(delimiters_mask), data.size-1)
        # super().__init__(data, new_lines, delimiters, header_data)
        # starts, ends = self._calculate_col_starts_and_ends(data, delimiters)
        # n_fields = next(i for i, d in enumerate(ends) if data[d] == '\n') + 1
        # self._n_cols = n_fields
        # self._wig_col_starts = starts.reshape(-1, n_fields)
        # self._wig_col_ends = ends.reshape(-1, n_fields)

    @classmethod
    def _get_buffer_extractor(cls, data, new_lines) -> TextBufferExtractor:
        delimiters_mask = (data == cls.DELIMITER)
        delimiters_mask[new_lines] = True
        delimiters = np.append(np.flatnonzero(delimiters_mask), data.size - 1)
        starts, ends = cls._calculate_col_starts_and_ends(data, delimiters)
        n_fields = next(i for i, d in enumerate(ends) if data[d] == '\n') + 1
        return TextBufferExtractor(data,
                                   starts.reshape(-1, n_fields),
                                   ends.reshape(-1, n_fields))

    @classmethod
    def from_raw_buffer(cls, chunk: np.ndarray, header_data=None) -> "DelimitedBuffer":
        """Make EncodedArray of the chunk and extract all complete lines

        Also find all delimiters in the buffer

        Parameters
        ----------
        chunk : np.ndarray
            Raw bytes chunk as array
        header_data : 6
            Any header data that was read in `read_header`

        Returns
        -------
        DelimitedBuffer
            DelimitedBuffer object of all complete lines

        """
        chunk = EncodedArray(chunk, BaseEncoding)
        new_lines = np.flatnonzero(chunk == '\n')
        extractor = cls._get_buffer_extractor(
            chunk[:new_lines[-1] + 1], new_lines[:-1])
        return cls(extractor, header_data)

    def __get_integers(self, cols: list) -> np.ndarray:
        """Get integers from integer string

        Extract integers from the specified columns

        Parameters
        ----------
        cols : list
            list of columns containing integers

        """
        assert np.all(cols < self._n_cols), (str(self._data), cols, self._n_cols)
        cols = np.asanyarray(cols)
        integer_starts = self._col_starts(cols)
        integer_ends = self._col_ends(cols)
        array = self._move_intervals_to_2d_array(integer_starts, integer_ends, fill_value='0')
        try:
            digits = as_encoded_array(array, DigitEncoding).raw()
        except EncodingError as e:
            row_number = e.offset // array.shape[-1]  # rows._shape.starts, e.offset, side="right")-1
            raise FormatException(e.args[0], line_number=row_number)
        powers = 10 ** np.arange(digits.shape[-1])[::-1]
        return digits.dot(powers).reshape(-1, cols.size)


class GFFBuffer(DelimitedBufferWithInernalComments):
    dataclass = GFFEntry
