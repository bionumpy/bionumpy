import sys
import numpy as np
from io import FileIO
from npstructures import RaggedArray, RaggedView, RaggedShape
from .exceptions import FormatException
from ..bnpdataclass import bnpdataclass
from ..datatypes import SequenceEntry, SequenceEntryWithQuality
from ..encoded_array import EncodedArray, EncodedRaggedArray
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
        self.size = self._data.size

    def raise_if(condition, *args, **kwargs):
        if condition:
            raise FormatException(*args, **kwargs)

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
        comment = cls.COMMENT
        if isinstance(comment, str):
            comment = ord(comment)
        for line in file_object:
            if line[0] != comment:
                file_object.seek(-len(line), 1)
                break

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
            # indices, shape = RaggedView(starts, lens).get_flat_indices()
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
        n_new_lines = sum(np.sum(EncodedArray(chunk, BaseEncoding) == NEWLINE) for chunk in chunks)
        return n_new_lines >= cls.n_lines_per_entry


class OneLineBuffer(FileBuffer):
    """ Base class for file formats where data fields are contained in lines."""

    n_lines_per_entry = 2
    _buffer_divisor = 32

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
        assert n_lines >= cls.n_lines_per_entry, "No complete entry in buffer. Try increasing chunk_size."
        new_lines = new_lines[: n_lines - (n_lines % cls.n_lines_per_entry)]
        chunk = chunk[: new_lines[-1] + 1]
        return cls(chunk[: new_lines[-1] + 1], new_lines)

    def get_data(self) -> bnpdataclass:
        """Get and parse fields from each line"""
        self.validate_if_not()
        starts = np.insert(self._new_lines, 0, -1)
        lengths = np.diff(starts)
        self.lines = EncodedRaggedArray(self._data, RaggedShape(lengths))
        sequences = self.lines[1::self.n_lines_per_entry, :-1]
        headers = self.lines[:: self.n_lines_per_entry, 1:-1]
        return SequenceEntry(headers, sequences)

    def count_entries(self) -> int:
        """Count number of entries in file"""
        return len(self._new_lines)//self.n_lines_per_entry

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

        name_lengths = entries.name.lengths
        sequence_lengths = entries.sequence.lengths
        line_lengths = np.hstack(
            (name_lengths[:, None] + 2, sequence_lengths[:, None] + 1)
        ).ravel()
        buf = EncodedArray(np.empty(line_lengths.sum(), dtype=np.uint8), BaseEncoding)
        lines = EncodedRaggedArray(buf, line_lengths)
        step = cls.n_lines_per_entry
        lines[0::step, 1:-1] = entries.name
        lines[1::step, :-1] = EncodedRaggedArray(
            EncodedArray(entries.sequence.encoding.decode(entries.sequence.ravel()),  entries.sequence.encoding),entries.sequence.shape)

        lines[0::step, 0] = ">"

        lines[:, -1] = "\n"
        return buf

    def _validate(self):
        n_lines = self._new_lines.size
        assert n_lines % self.n_lines_per_entry == 0, "Wrong number of lines in buffer"
        header_idxs = (
            self._new_lines[self.n_lines_per_entry - 1 : -1 : self.n_lines_per_entry]
            + 1
        )
        if np.any(self._data[header_idxs] != self.HEADER) or self._data[0] != self.HEADER:
            if self._data[0] != self.HEADER:
                line_number = 0
            else:
                line_number = (np.flatnonzero(self._data[header_idxs] != self.HEADER)[0]+1)*self.n_lines_per_entry
            raise FormatException(f"Expected header line to start with {self.HEADER}" % self._data, line_number=line_number)
        self._is_validated = True


class TwoLineFastaBuffer(OneLineBuffer):
    HEADER = ">"# 62
    n_lines_per_entry = 2
    dataclass = SequenceEntry


class FastQBuffer(OneLineBuffer):
    HEADER = "@"
    n_lines_per_entry = 4
    dataclass = SequenceEntryWithQuality

    def get_data(self):
        seq_entry = super().get_data()
        quality = self.lines[3 :: self.n_lines_per_entry, :-1]
        return SequenceEntryWithQuality(
            seq_entry.name, seq_entry.sequence,
            #quality)
            QualityEncoding.encode(quality))

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

    def _validate(self):
        super()._validate()
        if np.any(self._data[self._new_lines[1::self.n_lines_per_entry] + 1] != "+"):
            entry_number = np.flatnonzero(self._data[self._new_lines[1::self.n_lines_per_entry] + 1] != "+")[0]
            line_number = 2+entry_number*self.n_lines_per_entry
            raise FormatException(f"Expected '+' at third line of entry in {self._data}", line_number=line_number)
        
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
