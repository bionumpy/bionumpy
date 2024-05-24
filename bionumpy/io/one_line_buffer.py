import dataclasses
from typing import List, Tuple, Optional, Union

import numpy as np

from .dump_csv import get_column
from ..encoded_array import EncodedArray, BaseEncoding, change_encoding, EncodedRaggedArray
from ..bnpdataclass import bnpdataclass, BNPDataClass
from ..datatypes import SequenceEntry
from .exceptions import FormatException
from .file_buffers import FileBuffer, TextBufferExtractor, IncompleteEntryException, NEWLINE, TextThroughputExtractor


class OneLineBuffer(FileBuffer):
    """ Base class for file formats where data fields are contained in lines."""

    n_lines_per_entry = 2
    _buffer_divisor = 32
    _line_offsets = (1, 0)
    _empty_lines = []

    def __init__(self, buffer_extractor: TextBufferExtractor):
        self._buffer_extractor = buffer_extractor
        self._is_validated = True

    @property
    def n_lines(self) -> int:
        return len(self._buffer_extractor)*self.n_lines_per_entry

    @property
    def data(self) -> EncodedArray:
        return self._buffer_extractor.data


    @classmethod
    def contains_complete_entry(cls, chunks: List[EncodedArray]) -> Tuple[bool, EncodedArray]:
        if len(chunks) == 1:
            try:
                return True, cls.from_raw_buffer(chunks[0])
            except IncompleteEntryException:
                return False
        return super().contains_complete_entry(chunks)

    @classmethod
    def from_raw_buffer(cls, chunk: EncodedArray, header_data=None) -> "OneLineBuffer":
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

    def get_data(self) -> BNPDataClass:
        """Get and parse fields from each line"""
        headers, sequences = [self._buffer_extractor.get_field_by_number(i) for i in (0, 1)]
        return SequenceEntry(headers, sequences)

    def get_field_by_number(self, i: int, t: Optional[type]= None) -> Union[np.ndarray, EncodedArray, EncodedRaggedArray]:
        """ Get a field indexed by number"""

        text = self._buffer_extractor.get_field_by_number(i)
        if t is not None:
            return self._get_parser(t)(text)
        return text

    def __getitem__(self, idx):
        return self.__class__(self._buffer_extractor[idx])

    def count_entries(self) -> int:
        """Count number of entries in file"""
        return len(self._buffer_extractor)

    def get_field_range_as_text(self, start: int, end: int) -> EncodedRaggedArray:
        """Get a range of fields as text"""
        assert end == start+1
        return self.get_text_field_by_number(start)

    @classmethod
    def from_data(cls, entries: BNPDataClass) -> "OneLineBuffer":
        """Convert the data from the entries into a buffer that can be written to file

        Extract the name and the sequence and make a buffer with alternating lines

        Parameters
        ----------
        entries : BNPDataClass
            The entries to be written to the buffer


        Returns
        -------
        "OneLineBuffer"
            A ASCII encoded buffer
        """
        data_dict = [(field.type, getattr(entries, field.name)) for field in dataclasses.fields(entries)]
        columns= [get_column(value, key) for key, value in data_dict]
        return cls.join_fields(columns)

    @classmethod
    def join_fields(cls, fields: List[EncodedRaggedArray]) -> EncodedArray:
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

    def get_text_field_by_number(self, i: int) -> Union[np.ndarray, EncodedArray, EncodedRaggedArray]:
        return self.get_field_by_number(i)

    @classmethod
    def _get_buffer_extractor(cls, data, new_lines):
        tmp = np.insert(new_lines, 0, -1)+1
        field_ends = new_lines.reshape(-1, cls.n_lines_per_entry)
        field_ends = cls._modify_for_carriage_return(field_ends, data)
        field_starts = tmp[:-1].reshape(-1, cls.n_lines_per_entry)+(np.array(cls._line_offsets))
        entry_starts = tmp[:-1:cls.n_lines_per_entry]
        entry_ends = tmp[::cls.n_lines_per_entry][1:]

        return TextThroughputExtractor(data,
                                       field_starts,
                                       field_ends=field_ends,
                                       entry_starts=entry_starts,
                                       entry_ends=entry_ends)


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

    @classmethod
    def _modify_for_carriage_return(cls, field_ends, data):
        if field_ends.size == 0 or field_ends[0, 0] < 1:
            return field_ends
        last_chars = data[field_ends[:cls.n_lines_per_entry, 0] - 1]
        if not np.any(last_chars == "\r"):
            return field_ends
        return field_ends - (data[field_ends-1] == '\r')


class TwoLineFastaBuffer(OneLineBuffer):
    """
    Buffer for fasta files where each entry is contained in two lines (one for header and one for sequence)
    For multi-line fasta files, use MultiLineFastaBuffer
    """
    HEADER = ">"  # 62
    n_lines_per_entry = 2
    dataclass = SequenceEntry
