import numpy as np
from npstructures import RaggedArray, RaggedView, RaggedShape, npdataclass
from .encodings import BaseEncoding, QualityEncoding
from .sequences import to_ascii, ASCIIText, Sequences
from .datatypes import SequenceEntry, SequenceEntryWithQuality

NEWLINE = 10


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

    def __init__(self, data, new_lines):
        self._data = np.asanyarray(data).view(ASCIIText)
        self._new_lines = np.asanyarray(new_lines)
        self._is_validated = False
        self.size = self._data.size

    @classmethod
    def read_header(cls, file_object):
        pass

    @classmethod
    def from_raw_buffer(cls, raw_buffer, header_data=None) -> "FileBuffer":
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
    def from_data(cls, data: npdataclass) -> "FileBuffer":
        """Create FileBuffer from a data set

        Create a FileBuffer that can be written to file

        Parameters
        ----------
        data : npdataclass
            Data set containing the data to be written

        Returns
        -------
        'FileBuffer'
            FileBuffer containing the data
        """
        return NotImplemented

    def validate_if_not(self):
        if not self._is_validated:
            self._validate()

    def get_data(self) -> npdataclass:
        """Extract the data from the buffer

        The default way to extract data from the the buffer

        Returns
        -------
        npdataclass
            Data set containing the data from the buffer
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
        indices, shape = RaggedView(starts, lens).get_flat_indices()
        return Sequences(self._data[indices], shape)

    def _move_2d_array_to_intervals(self, array, starts, ends):
        n_chars = ends - starts
        n_intervals = starts.size
        max_chars = array.shape[-1]
        to_indices = ends[::-1, None]-max_chars+np.arange(max_chars)
        self._data[to_indices] = array[::-1]


class OneLineBuffer(FileBuffer):
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
        new_lines = np.flatnonzero(chunk == NEWLINE)
        n_lines = new_lines.size
        assert n_lines >= cls.n_lines_per_entry, "No complete entry in buffer. Try increasing chunk_size."
        new_lines = new_lines[: n_lines - (n_lines % cls.n_lines_per_entry)]
        return cls(chunk[: new_lines[-1] + 1], new_lines)

    def get_sequences(self) -> RaggedArray:
        self.validate_if_not()
        sequence_starts = self._new_lines[:: self.n_lines_per_entry] + 1
        sequence_lens = self._new_lines[1 :: self.n_lines_per_entry] - sequence_starts
        indices, shape = RaggedView(sequence_starts, sequence_lens).get_flat_indices()
        m = indices.size
        d = m % self._buffer_divisor
        seq = np.empty(m - d + self._buffer_divisor, dtype=self._data.dtype).view(ASCIIText)
        seq[:m] = self._data[indices]
        return RaggedArray(seq, shape)

    def get_data(self):
        self.validate_if_not()
        starts = np.insert(self._new_lines, 0, -1)
        lengths = np.diff(starts)
        self.lines = RaggedArray(self._data, RaggedShape(lengths))
        sequences = self.lines[1 :: self.n_lines_per_entry, :-1]
        headers = self.lines[:: self.n_lines_per_entry, 1:-1]
        return SequenceEntry(headers, sequences)

    def count_entries(self):
        return len(self._new_lines)//self.n_lines_per_entry

    @classmethod
    def from_data(cls, entries):
        name_lengths = entries.name.shape.lengths
        sequence_lengths = entries.sequence.shape.lengths
        line_lengths = np.hstack(
            (name_lengths[:, None] + 2, sequence_lengths[:, None] + 1)
        ).ravel()
        buf = np.empty(line_lengths.sum(), dtype=np.uint8).view(ASCIIText)
        lines = RaggedArray(buf, line_lengths)
        step = cls.n_lines_per_entry
        lines[0::step, 1:-1] = entries.name
        lines[1::step, :-1] = entries.sequence

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
        assert np.all(self._data[header_idxs] == self.HEADER)
        self._is_validated = True


class TwoLineFastaBuffer(OneLineBuffer):
    HEADER = 62
    n_lines_per_entry = 2
    _encoding = BaseEncoding
    dataclass = SequenceEntry


class FastQBuffer(OneLineBuffer):
    HEADER = "@"
    n_lines_per_entry = 4
    dataclass = SequenceEntryWithQuality

    def get_data(self):
        seq_entry = super().get_data()
        quality = self.lines[3 :: self.n_lines_per_entry, :-1]
        return SequenceEntryWithQuality(
            seq_entry.name, seq_entry.sequence, quality)

    @classmethod
    def _get_line_lens(cls, entries):
        name_lengths = entries.name.shape.lengths[:, None]
        sequence_lengths = entries.sequence.shape.lengths[:, None]
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
    def from_data(cls, entries):
        line_lengths = cls._get_line_lens(entries)
        buf = np.empty(line_lengths.sum(), dtype=np.uint8).view(ASCIIText)
        lines = RaggedArray(buf, line_lengths)
        step = cls.n_lines_per_entry
        lines[0::step, 1:-1] = entries.name
        lines[1::step, :-1] = entries.sequence
        lines[2::step, 0] = "+"
        lines[3::step, :-1] = to_ascii(entries.quality, QualityEncoding)
        lines[0::step, 0] = cls.HEADER
        lines[:, -1] = "\n"

        return buf
