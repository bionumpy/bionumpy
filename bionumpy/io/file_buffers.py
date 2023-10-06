import sys
from typing import List

import numpy as np
from io import FileIO
from npstructures import RaggedArray, RaggedView, RaggedShape
from .exceptions import FormatException
from ..bnpdataclass import bnpdataclass
from ..datatypes import SequenceEntry, SequenceEntryWithQuality
from ..encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from ..encodings import QualityEncoding, BaseEncoding

NEWLINE = "\n"


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
    def size(self):
        return self.data.size

    @property
    def data(self):
        return self._data

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
        return len(self._new_lines)

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
        max_chars = np.max(ends-starts)
        view_starts = (ends-max_chars)
        indices = view_starts[:, None]+np.arange(max_chars)
        array = self._data[indices.ravel()]
        zeroed, _ = RaggedView(np.arange(starts.size)*max_chars, max_chars-(ends-starts)).get_flat_indices()
        array[zeroed] = fill_value
        return array.reshape((-1, max_chars))

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

    #def get_field_by_number(self, field_nr: int, field_type: type=object):
    #    raise NotImplementedError


class IncompleteEntryException(Exception):
    pass


class TextBufferExtractor:
    def __init__(self, data: EncodedArray, field_starts: np.ndarray, field_ends: np.ndarray):
        '''
        field_starts: n_entries x n_fields
        field_ends: n_entries x n_fields
        '''
        self._data = data
        self._field_starts = field_starts
        self._field_ends = field_ends
        self._field_lens = field_ends-field_starts
        self._n_fields = field_starts.shape[1]

    @property
    def data(self):
        return self._data

    def __len__(self):
        return len(self._field_starts)

    def __getitem__(self, idx):
        return self.__class__(self._data, field_starts=self._field_starts[idx], field_ends=self._field_ends[idx])

    def get_field_by_number(self, field_nr: int):
        e = EncodedRaggedArray(self._data, RaggedView(self._field_starts.ravel(), self._field_lens.ravel()))
        return e[field_nr::self._n_fields]


class OneLineBuffer(FileBuffer):
    """ Base class for file formats where data fields are contained in lines."""

    n_lines_per_entry = 2
    _buffer_divisor = 32
    _line_offsets = (1, 0)
    _empty_lines = []

    def __init__(self, buffer_extractor: TextBufferExtractor):
        # super().__init__(*args, **kwargs)
        self._is_validated=False
        self._buffer_extractor = buffer_extractor
        self._is_validated = True

    @property
    def n_lines(self):
        return len(self._buffer_extractor)*self.n_lines_per_entry

    @property
    def data(self):
        return self._buffer_extractor.data

    @classmethod
    def _get_buffer_extractor(cls, data, new_lines):
        tmp = np.insert(new_lines, 0, -1)
        ends = tmp[1:].reshape(-1, cls.n_lines_per_entry)
        starts = tmp[:-1].reshape(-1, cls.n_lines_per_entry)+(np.array(cls._line_offsets)+1)
        return TextBufferExtractor(data, starts, ends)

    @classmethod
    def contains_complete_entry(cls, chunks):
        if len(chunks) == 1:
            try:
                return True, cls.from_raw_buffer(chunks[0])
            except IncompleteEntryException:
                return False
        return super().contains_complete_entry(chunks)

    @classmethod
    def from_raw_buffer(cls, chunk, header_data=None) -> "OneLineBuffer":
        """Create a buffer with full entries

        Extract complete entries, i. e. a number of lines that is divisible by lines per entry

        Parameters
        ----------
        chunk : np.ndarray
            Raw buffer with data that might end with incomplete entry

        Returns
        -------
        'OneLineBuffer'
            Buffer with complete entries

        Examples
        --------
        8

        """
        assert header_data is None
        chunk = EncodedArray(chunk, BaseEncoding)
        new_lines = np.flatnonzero(chunk == NEWLINE)
        n_lines = new_lines.size
        if n_lines < cls.n_lines_per_entry:
            raise IncompleteEntryException("No complete entry in buffer. Try increasing chunk_size.")
        new_lines = new_lines[: n_lines - (n_lines % cls.n_lines_per_entry)]
        chunk = chunk[: new_lines[-1] + 1]
        data = chunk[: new_lines[-1] + 1]
        cls._validate(data, new_lines)
        return cls(cls._get_buffer_extractor(data, new_lines))

    @property
    def lines(self):
        if not hasattr(self, "__lines"):
            starts = np.insert(self._new_lines, 0, -1)
            lengths = np.diff(starts)
            self.__lines = EncodedRaggedArray(self._data, RaggedShape(lengths))
        return self.__lines

    @property
    def entries(self):
        if not hasattr(self, "__entries"):
            lengths = np.diff(self._new_lines[self.n_lines_per_entry-1::self.n_lines_per_entry])
            lengths = np.insert(lengths, 0, self._new_lines[self.n_lines_per_entry-1]+1)
            self.__entries = EncodedRaggedArray(self._data, RaggedShape(lengths))
        return self.__entries

    def get_data(self) -> bnpdataclass:
        """Get and parse fields from each line"""
        headers, sequences = [self._buffer_extractor.get_field_by_number(i) for i in (0, 1)]
        return SequenceEntry(headers, sequences)

    def get_field_by_number(self, i: int, t: type=object):
        """ Get a field indexed by number"""
        return self._buffer_extractor.get_field_by_number(i)

    def __getitem__(self, idx):
        return self.__class__(self._buffer_extractor[idx])

    def count_entries(self) -> int:
        """Count number of entries in file"""
        return len(self._buffer_extractor)

    def get_field_range_as_text(self, start, end):
        """Get a range of fields as text"""
        assert end == start+1
        return self.get_text_field_by_number(start)

    @classmethod
    def from_data(cls, entries: bnpdataclass) -> "OneLineBuffer":
        """Convert the data from the entries into a buffer that can be written to file
        
        Extract the name and the sequence and make a buffer with alternating lines

        Parameters
        ----------
        entries : bnpdataclass
            The entries to be written to the buffer
        

        Returns
        -------
        "OneLineBuffer"
            A ASCII encoded buffer
        """

        names = entries.name
        sequences = entries.sequence
        name_lengths = names.lengths
        sequence_lengths = sequences.lengths
        line_lengths = np.hstack(
            (name_lengths[:, None] + 2, sequence_lengths[:, None] + 1)
        ).ravel()
        buf = EncodedArray(np.empty(line_lengths.sum(), dtype=np.uint8), BaseEncoding)
        lines = EncodedRaggedArray(buf, line_lengths)
        step = cls.n_lines_per_entry
        lines[0::step, 1:-1] = names
        lines[1::step, :-1] = EncodedRaggedArray(
            EncodedArray(sequences.encoding.decode(sequences.ravel()), sequences.encoding), sequences.shape)

        lines[0::step, 0] = ">"

        lines[:, -1] = "\n"
        return buf

    @classmethod
    def join_fields(cls, fields: List[EncodedRaggedArray]):
        field_lengths = np.hstack([field.shape[1][:, None] for field in fields])
        line_lengths = field_lengths+1
        for i in range(len(fields)):
            line_lengths[:, i] += cls._line_offsets[i]
        entry_lengths = line_lengths.sum(axis=-1)
        buffer_size = entry_lengths.sum()
        buf = EncodedArray(np.empty(buffer_size, dtype=np.uint8), BaseEncoding)
        lines = EncodedRaggedArray(buf, line_lengths.ravel())
        step = cls.n_lines_per_entry
        for i, field in enumerate(fields):
            lines[i::step, cls._line_offsets[i]:-1] = field
        lines[0::step, 0] = cls.HEADER
        lines[:, -1] = "\n"
        return buf

    @classmethod
    def _validate(cls, data, new_lines):
        header = cls.HEADER
        if data.size == 0 and new_lines.size==0:
            return
        n_lines = new_lines.size
        n_lines_per_entry = cls.n_lines_per_entry
        assert n_lines % n_lines_per_entry == 0, "Wrong number of lines in buffer"
        header_idxs = (
                new_lines[n_lines_per_entry - 1: -1: n_lines_per_entry]
                + 1
        )

        if np.any(data[header_idxs] != header) or data[0] != header:
            if data[0] != header:
                line_number = 0
            else:
                line_number = (np.flatnonzero(data[header_idxs] != header)[0] + 1) * n_lines_per_entry
            raise FormatException(f"Expected header line to start with {header}" % data, line_number=line_number)

    def get_text_field_by_number(self, i):
        return self.get_field_by_number(i)


class TwoLineFastaBuffer(OneLineBuffer):
    HEADER = ">"# 62
    n_lines_per_entry = 2
    dataclass = SequenceEntry


class FastQBuffer(OneLineBuffer):
    HEADER = "@"
    n_lines_per_entry = 4
    dataclass = SequenceEntryWithQuality
    _line_offsets = (1, 0, 0, 0)
    _empty_lines = [2]

    def get_text_field_by_number(self, i: int) -> EncodedRaggedArray:
        if i == 2:
            return self._buffer_extractor.get_field_by_number(3)
        return super().get_text_field_by_number(i)

    def get_field_by_number(self, i: int, t: type=object):
        if i == 2:
            return QualityEncoding.encode(self.get_text_field_by_number(i))
        else:
            return super().get_field_by_number(i, t)

    def get_data(self):
        seq_entry = super().get_data()
        quality = self.get_field_by_number(2, QualityEncoding)
        return SequenceEntryWithQuality(
            seq_entry.name, seq_entry.sequence,quality)

    @classmethod
    def _get_line_lens(cls, entries):
        name_lengths = entries.name.lengths[:, None]
        sequence_lengths = entries.sequence.lengths[:, None]
        return (
            np.hstack(
                (
                    name_lengths + 1,
                    sequence_lengths,
                    np.ones_like(sequence_lengths),
                    sequence_lengths,
                )
            ).ravel()
            + 1
        )

    @classmethod
    def _validate(cls, data, new_lines):
        super()._validate(data, new_lines)
        n_lines_per_entry = cls.n_lines_per_entry
        if np.any(data[new_lines[1::n_lines_per_entry] + 1] != "+"):
            entry_number = np.flatnonzero(data[new_lines[1::n_lines_per_entry] + 1] != "+")[0]
            line_number = 2 + entry_number * n_lines_per_entry
            raise FormatException(f"Expected '+' at third line of entry in {data}", line_number=line_number)

    @classmethod
    def join_fields(cls, fields: List[EncodedRaggedArray]):
        plus_line = as_encoded_array(['+']*len(fields[0]))
        return super().join_fields(fields[:2]+[plus_line]+fields[2:])

    @classmethod
    def from_data(cls, entries):
        line_lengths = cls._get_line_lens(entries)
        buf = EncodedArray(np.empty(line_lengths.sum(), dtype=np.uint8), BaseEncoding)
        lines = EncodedRaggedArray(buf, line_lengths)
        step = cls.n_lines_per_entry
        lines[0::step, 1:-1] = entries.name
        lines[1::step, :-1] = EncodedRaggedArray(EncodedArray(entries.sequence.encoding.decode(entries.sequence.ravel()), BaseEncoding), entries.sequence.shape)
        lines[2::step, 0] = "+"
        lines[3::step, :-1] = EncodedArray(QualityEncoding.decode(entries.quality.ravel()), BaseEncoding)
        lines[0::step, 0] = cls.HEADER
        lines[:, -1] = "\n"

        return buf
