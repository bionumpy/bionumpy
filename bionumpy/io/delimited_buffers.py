import io
import logging
import dataclasses
from typing import List
from npstructures import RaggedArray, ragged_slice
from ..bnpdataclass import bnpdataclass, BNPDataClass
from ..datatypes import (Interval, VCFGenotypeEntry,
                         SequenceEntry, VCFEntry, Bed12, Bed6,
                         GFFEntry, SAMEntry, ChromosomeSize, NarrowPeak, PhasedVCFGenotypeEntry,
                         GfaPath)
from ..encoded_array import EncodedArray, EncodedRaggedArray
from ..encoded_array import as_encoded_array
from ..encodings import Encoding
from ..encodings.exceptions import EncodingError
from ..encodings.vcf_encoding import GenotypeRowEncoding, PhasedGenotypeRowEncoding
from ..encodings.alphabet_encoding import get_alphabet_encodings, DigitEncoding
from ..encoded_array import get_base_encodings, BaseEncoding
from ..util import is_subclass_or_instance
from .file_buffers import FileBuffer, NEWLINE
from .strops import (
    ints_to_strings, split, str_to_int, str_to_float,
    int_lists_to_strings, float_to_strings)
from .dump_csv import dump_csv
from .exceptions import FormatException
import numpy as np


class DelimitedBuffer(FileBuffer):
    """Base class for file buffers for delimited files such as csv or tsv.
    Each line should correspond to an entry, and each column to a variable.

    Provides convenience methods for extracting and decoding integers from columns,
    and text from columns into Sequences objects
    """

    DELIMITER = "\t"
    COMMENT = "#"
    INCLUDE_HEADER = False
    n_lines_per_entry = 1

    def __init__(self, data: EncodedArray, new_lines: np.ndarray, delimiters: np.ndarray = None, header_data=None):
        super().__init__(data, new_lines)
        if delimiters is None:
            delimiters = np.concatenate(
                ([-1], np.flatnonzero(self._data == self.DELIMITER), self._new_lines)
            )
            delimiters.sort(kind="mergesort")
        self._delimiters = delimiters
        self._header_data = header_data

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
        n_fields = next((i + 1 for i, v in enumerate(delimiters) if chunk[v] == "\n"), None)
        if n_fields is None:
            logging.warning("Foud no new lines. Chunk size may be too low. Try increasing")
            raise 
        new_lines = delimiters[(n_fields - 1)::n_fields]
        #if not all(chunk[new_lines] == "\n"):
        #     raise FormatException("Irregular number of columns in file")
        delimiters = np.concatenate(([-1], delimiters[:n_fields * len(new_lines)]))
        return cls(chunk[:new_lines[-1] + 1], new_lines, delimiters, header_data)

    def get_integers(self, cols: list) -> np.ndarray:
        """Get integers from integer string

        Extract integers from the specified columns

        Parameters
        ----------
        cols : list
            list of columns containing integers

        """
        assert np.all(cols < self._n_cols), (str(self._data), cols, self._n_cols)
        cols = np.asanyarray(cols)
        integer_starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, cols] + 1
        integer_ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, cols]
        integers = self._extract_integers(integer_starts.ravel(), integer_ends.ravel())
        return integers.reshape(-1, cols.size)

    def get_floats(self, cols) -> np.ndarray:
        """Get floats from float string

        Extract floats from the specified columns

        Parameters
        ----------
        cols : list
            list of columns containing integers

        Examples
        --------
        FIXME: Add docs.

        """
        cols = np.asanyarray(cols)
        float_starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, cols] + 1
        float_ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, cols]
        floats = str_to_float(ragged_slice(self._data, float_starts, float_ends))
        return floats.reshape(-1, cols.size)

    def get_text(self, col, fixed_length=True, keep_sep=False):
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
        self.validate_if_not()
        assert np.max(col) < self._n_cols, (col, self._n_cols, self._data, self.__class__, self.dataclass)
        starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, col] + 1
        ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, col]
        if keep_sep:
            ends += 1
        if fixed_length:
            return self._move_intervals_to_2d_array(starts, ends)
        else:
            return self._move_intervals_to_ragged_array(starts, ends)

    def get_column_range_as_text(self, col_start, col_end, keep_sep=False):
        """Get multiple columns as text

       Parameters
        ----------
        col_start : int
            column index to start at
        col_end : int
            column index to end at
        keep_sep : bool
            keep seperator at end
        """
        self.validate_if_not()
        assert col_start < col_end
        assert col_start < self._n_cols, self._n_cols
        assert col_end <= self._n_cols

        starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, col_start] + 1
        ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, col_end-1]
        if keep_sep:
            ends += 1

        return self._move_intervals_to_ragged_array(starts, ends)

    def get_text_range(self, col, start=0, end=None) -> np.ndarray:
        """Get substrings of a column

        Extract the text from start to end of each entry in column

        Parameters
        ----------
        col : int
            column index
        start : int
            start of substring
        end : int
            end of substring

        Returns
        -------
        np.ndarray
            array containing the extracted substrings

        Examples
        --------
        FIXME: Add docs.

        """
        self.validate_if_not()
        starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, col] + 1 + start
        if end is not None:
            return self._data[starts[..., np.newaxis] + np.arange(end - start)].reshape(-1, end - start)
        else:
            ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, col]
        return self._move_intervals_to_2d_array(starts.ravel(), ends.ravel())

    def _extract_integers(self, integer_starts, integer_ends):
        rows = self._data[integer_starts:integer_ends]
        try:
            strs = str_to_int(rows)
        except EncodingError as e:
            row_number = np.searchsorted(rows._shape.starts, e.offset, side="right")-1
            raise FormatException(e.args[0], line_number=row_number)
        return strs

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
        should_be_new_lines = chunk[delimiters[n_delimiters_per_line-1::n_delimiters_per_line]]
        if delimiters.size % n_delimiters_per_line != 0 or np.any(should_be_new_lines != "\n"):
            offending_line = np.flatnonzero(should_be_new_lines != "\n")[0]
            lines = split(self._data, '\n')
            raise FormatException(f"Irregular number of delimiters per line ({delimiters.size}, {n_delimiters_per_line}): {lines}", line_number=offending_line)
        self._validated = True

    @classmethod
    def from_data(cls, data: BNPDataClass) -> "DelimitedBuffer":
        """Put each field of the dataclass into a column in a buffer.

        Parameters
        data : bnpdataclass
            Data
        """
        data_dict = [(field.type, getattr(data, field.name)) for field in dataclasses.fields(data)]
        return dump_csv(data_dict, cls.DELIMITER)

    @classmethod
    def make_header(cls, data: bnpdataclass):
        """makes a header from field names separated by delimiter"""
        return bytes(cls.DELIMITER.join([field.name for field in dataclasses.fields(data)]) + "\n", 'ascii')

    def get_data(self) -> BNPDataClass:
        """Parse the data in the buffer according to the fields in _dataclass

        Returns
        -------
        _dataclass
            Dataclass with parsed data

        """
        self.validate_if_not()
        columns = {}
        fields = dataclasses.fields(self.dataclass)
        for col_number, field in enumerate(fields):
            if field.type is None:
                col = None
            elif field.type == str or is_subclass_or_instance(field.type, Encoding):
                col = self.get_text(col_number, fixed_length=False)
            elif field.type == int:
                col = self.get_integers(col_number).ravel()
            elif field.type == float:
                col = self.get_floats(col_number).ravel()
            elif field.type == -1:
                col = self.get_integers(col_number).ravel() - 1
            elif field.type == List[int]:
                col = self.get_split_ints(col_number)
            elif field.type == List[bool]:
                col = self.get_split_ints(col_number, sep="").astype(bool)
            else:
                assert False, field
            columns[field.name] = col
        n_entries = len(next(col for col in columns if col is not None))
        columns = {c: value if c is not None else np.empty((n_entries, 0))
                   for c, value in columns.items()}
        return self.dataclass(**columns)

    def get_split_ints(self, col: int, sep: str = ",") -> RaggedArray:
        """Split a column of separated integers into a raggedarray

        Parameters
        ----------
        col : int
            Column number
        sep : 5
            Characted used to separate the integers

        Examples
        --------
        7

        """
        self.validate_if_not()
        starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, col] + 1
        ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, col] + len(sep)
        text = self._data[starts:ends]
        if len(sep):
            text[:, -1] = sep
            int_strings = split(text.ravel()[:-1], sep=sep)
            return RaggedArray(str_to_int(int_strings), (text == sep).sum(axis=-1))
        else:
            mask = as_encoded_array(text.ravel(), DigitEncoding).raw()
            return RaggedArray(mask, text.shape)

    def count_entries(self) -> int:
        """Count the number of entries in the buffer"""
        return len(self._new_lines)


class GfaSequenceBuffer(DelimitedBuffer):
    dataclass = SequenceEntry

    def get_data(self):
        ids = self.get_text(1, fixed_length=False)
        sequences = self.get_text(col=2, fixed_length=False)
        return SequenceEntry(ids, sequences)

    @classmethod
    def from_data(cls, data: SequenceEntry) -> EncodedArray:
        return dump_csv([(str, as_encoded_array(["S"]*len(data))),
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
        directions = np.where(all_node_texts[:, -1]=="+", 1, -1)
        node_ids = RaggedArray(node_ids, lengths)
        directions = RaggedArray(directions, lengths)
        data =  GfaPath(name, node_ids, directions)
        return data



def get_bufferclass_for_datatype(_dataclass: bnpdataclass, delimiter: str = "\t", has_header: bool = False, comment: str = "#",
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
        INCLUDE_HEADER = has_header
        dataclass = _dataclass
        fields = None

        def __init__(self, data: EncodedArray, new_lines: np.ndarray, delimiters: np.ndarray = None, header_data: List[str] = None):
            super().__init__(data, new_lines, delimiters, header_data)
            self.set_fields_from_header(header_data)

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

        def set_fields_from_header(self, columns: List[str]):
            if not has_header:
                return None
            fields = dataclasses.fields(self.dataclass)
            self.fields = [next(field for field in fields if field.name == col) for col in columns]
            assert np.array_equal(columns, [field.name for field in self.fields])

        def get_data(self) -> _dataclass:
            """Parse the data in the buffer according to the fields in _dataclass

            Returns
            -------
            _dataclass
                Dataclass with parsed data

            """
            self.validate_if_not()
            columns = {}
            fields = self.fields if self.fields is not None else dataclasses.fields(self.dataclass)
            for col_number, field in enumerate(fields):
                if field.type is None:
                    col = None
                elif field.type == str or is_subclass_or_instance(field.type, Encoding):
                    col = self.get_text(col_number, fixed_length=False)
                elif field.type == int:
                    col = self.get_integers(col_number).ravel()
                elif field.type == float:
                    col = self.get_floats(col_number).ravel()
                elif field.type == List[int]:
                    col = self.get_split_ints(col_number)
                elif field.type == List[bool]:
                    col = self.get_split_ints(col_number, sep="").astype(bool)
                else:
                    assert False, field
                columns[field.name] = col
            n_entries = len(next(col for col in columns if col is not None))
            columns = {c: value if c is not None else np.empty((n_entries, 0))
                       for c, value in columns.items()}
            return self.dataclass(**columns)

    DatatypeBuffer.__name__ = _dataclass.__name__ + "Buffer"
    DatatypeBuffer.__qualname__ = _dataclass.__qualname__ + "Buffer"
    return DatatypeBuffer


class BedBuffer(DelimitedBuffer):
    dataclass = Interval

    def get_integers(self, cols: list) -> np.ndarray:
        """Get integers from integer string

        Extract integers from the specified columns

        Parameters
        ----------
        cols : list
            list of columns containing integers

        """
        assert np.all(cols < self._n_cols), (str(self._data), cols, self._n_cols)
        cols = np.asanyarray(cols)
        integer_starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, cols] + 1
        integer_ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, cols]
        array = self._move_intervals_to_2d_array(integer_starts, integer_ends, fill_value='0')
        try:
            digits = as_encoded_array(array, DigitEncoding).raw()
        except EncodingError as e:
            row_number = e.offset//array.shape[-1]#  rows._shape.starts, e.offset, side="right")-1
            raise FormatException(e.args[0], line_number=row_number)
        powers = 10**np.arange(digits.shape[-1])[::-1]
        return digits.dot(powers).reshape(-1, cols.size)
        #integers = self._extract_integers(integer_starts.ravel(), integer_ends.ravel())
        #return integers.reshape(-1, cols.size)


class Bed6Buffer(BedBuffer):
    dataclass = Bed6


class Bed12Buffer(Bed6Buffer):
    dataclass = Bed12



class VCFBuffer(DelimitedBuffer):
    dataclass = VCFEntry

    def get_data(self):
        data = super().get_data()
        data.position -= 1
        return data

    @classmethod
    def from_data(cls, data: BNPDataClass) -> "DelimitedBuffer":
        data = dataclasses.replace(data, position=data.position+1)
        return super().from_data(data)


class VCFMatrixBuffer(VCFBuffer):
    dataclass = VCFEntry
    genotype_dataclass = VCFGenotypeEntry
    genotype_encoding = GenotypeRowEncoding

    @classmethod
    def read_header(cls, file_object):
        prev_line = None
        for line in file_object:
            line = line.decode()
            if line[0] != cls.COMMENT:
                file_object.seek(-len(line), 1)
                break
            prev_line = line
        if prev_line is None:
            return []

        sample_names = prev_line.split("\t")[9:]
        return sample_names

    def get_data(self):
        data = super().get_data()
        genotypes = self.get_column_range_as_text(9, self._n_cols, keep_sep=True)
        genotypes = EncodedArray(self.genotype_encoding.encode(genotypes), self.genotype_encoding)
        return self.genotype_dataclass(*data.shallow_tuple(), genotypes)


class PhasedVCFMatrixBuffer(VCFMatrixBuffer):
    dataclass = VCFEntry
    genotype_dataclass = PhasedVCFGenotypeEntry
    genotype_encoding = PhasedGenotypeRowEncoding


class NarrowPeakBuffer(DelimitedBuffer):
    dataclass = NarrowPeak


class GFFBuffer(DelimitedBuffer):
    dataclass = GFFEntry


class SAMBuffer(DelimitedBuffer):
    dataclass = SAMEntry
    COMMENT = "@"


class ChromosomeSizeBuffer(DelimitedBuffer):
    dataclass = ChromosomeSize
