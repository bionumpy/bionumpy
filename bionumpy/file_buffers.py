import numpy as np
from npstructures import RaggedArray, RaggedView, RaggedShape, npdataclass
from .encodings import BaseEncoding, QualityEncoding
from .datatypes import SequenceEntry, SequenceEntryWithQuality
from .sequences import Sequences

NEWLINE = 10


class FileBuffer:
    _buffer_divisor = 1
    COMMENT = 0

    def __init__(self, data, new_lines):
        self._data = np.asanyarray(data)
        self._new_lines = np.asanyarray(new_lines)
        self._is_validated = False
        self.size = self._data.size

    @classmethod
    def from_raw_buffer(cls, raw_buffer) -> "FileBuffer":
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
        n_intervals = starts.size
        n_chars = ends - starts
        from_indices, _ = RaggedView(starts, n_chars).get_flat_indices()
        max_chars = np.max(n_chars)
        array = np.full(n_intervals * max_chars, fill_value, dtype=np.uint8)
        to_indices, _ = RaggedView(
            max_chars * np.arange(1, n_intervals + 1) - n_chars, n_chars
        ).get_flat_indices()
        array[to_indices] = self._data[from_indices]
        return array.reshape((n_intervals, max_chars))

    def _move_intervals_to_ragged_array(self, starts, ends=None, lens=None):
        if lens is None:
            lens = ends - starts
        indices, shape = RaggedView(starts, lens).get_flat_indices()
        return Sequences(self._data[indices], shape)


class OneLineBuffer(FileBuffer):
    n_lines_per_entry = 2
    _buffer_divisor = 32

    @classmethod
    def from_raw_buffer(cls, chunk) -> "OneLineBuffer":
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
        new_lines = np.flatnonzero(chunk == NEWLINE)
        n_lines = new_lines.size
        assert n_lines >= cls.n_lines_per_entry, "No complete entry in buffer"
        new_lines = new_lines[: n_lines - (n_lines % cls.n_lines_per_entry)]
        return cls(chunk[: new_lines[-1] + 1], new_lines)

    def get_sequences(self) -> Sequences:
        self.validate_if_not()
        sequence_starts = self._new_lines[:: self.n_lines_per_entry] + 1
        sequence_lens = self._new_lines[1 :: self.n_lines_per_entry] - sequence_starts
        indices, shape = RaggedView(sequence_starts, sequence_lens).get_flat_indices()
        m = indices.size
        d = m % self._buffer_divisor
        seq = np.empty(m - d + self._buffer_divisor, dtype=self._data.dtype)
        seq[:m] = self._data[indices]
        return Sequences(seq, shape)

    def get_data(self):
        self.validate_if_not()
        starts = np.insert(self._new_lines, 0, -1)
        lengths = np.diff(starts)
        self.lines = Sequences(self._data, RaggedShape(lengths))
        sequences = self.lines[1 :: self.n_lines_per_entry, :-1]
        headers = self.lines[:: self.n_lines_per_entry, 1:-1]
        return SequenceEntry(headers, sequences)

    @classmethod
    def from_data(cls, entries):
        name_lengths = entries.name.shape.lengths
        sequence_lengths = entries.sequence.shape.lengths
        line_lengths = np.hstack(
            (name_lengths[:, None] + 2, sequence_lengths[:, None] + 1)
        ).ravel()
        buf = np.empty(line_lengths.sum(), dtype=np.uint8)
        lines = RaggedArray(buf, line_lengths)
        step = cls.n_lines_per_entry
        lines[0::step, 1:-1] = entries.name
        lines[1::step, :-1] = entries.sequence

        lines[0::step, 0] = ord(">")

        lines[:, -1] = ord("\n")
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


class FastQBuffer(OneLineBuffer):
    HEADER = 64
    n_lines_per_entry = 4
    _encoding = BaseEncoding
    dataclass = SequenceEntryWithQuality

    def get_data(self):
        seq_entry = super().get_data()
        quality = QualityEncoding.encode(
            self.lines[3 :: self.n_lines_per_entry, :-1]
        )
        return SequenceEntryWithQuality(seq_entry.name, seq_entry.sequence, quality)

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
        buf = np.empty(line_lengths.sum(), dtype=np.uint8)
        lines = RaggedArray(buf, line_lengths)
        step = cls.n_lines_per_entry
        lines[0::step, 1:-1] = entries.name
        lines[1::step, :-1] = entries.sequence
        lines[2::step, 0] = ord("+")
        lines[3::step, :-1] = QualityEncoding.decode(entries.quality)
        lines[0::step, 0] = cls.HEADER
        lines[:, -1] = ord("\n")

        return buf
