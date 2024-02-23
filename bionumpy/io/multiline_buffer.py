from npstructures import RaggedArray
from .file_buffers import FileBuffer
from ..bnpdataclass import bnpdataclass
from ..encoded_array import EncodedArray, EncodedRaggedArray, encoded_array_from_nparray
from ..datatypes import SequenceEntry
import numpy as np

from ..encodings import BaseEncoding


class MultiLineBuffer(FileBuffer):
    SKIP_LAZY = True


class MultiLineFastaBuffer(MultiLineBuffer):
    _new_entry_marker = ">"
    n_characters_per_line = 80
    dataclass = SequenceEntry

    def __init__(self, data, new_lines, new_entries):
        super().__init__(data, new_lines)
        self._new_entries = new_entries

    @property
    def data(self):
        return self._data

    @property
    def n_lines(self):
        return len(self._new_lines)

    @classmethod
    def contains_complete_entry(cls, chunks):
        ends_with_new_line = False
        for chunk in chunks:
            chunk = EncodedArray(chunk, BaseEncoding)
            new_lines = np.flatnonzero(chunk[:-1] == "\n")
            new_entries = np.flatnonzero(chunk[new_lines+1] == cls._new_entry_marker)
            if new_entries.size >= 1:
                return True
            if ends_with_new_line and chunk[0] == cls._new_entry_marker:
                return True
            ends_with_new_line = chunk[-1] == "\n"
        return False

    def get_data(self):
        self.validate_if_not()
        line_starts = np.insert(self._new_lines + 1, 0, 0)
        line_ends = np.append(self._new_lines, self._data.size-1)
        line_ends = self._modify_ends_for_carriage_returns(line_ends, self._data)
        data = self._move_intervals_to_ragged_array(line_starts, line_ends)
        data.ravel()
        new_entries = np.insert(self._new_entries+1, 0, 0)
        n_lines_per_entry = np.diff(np.append(new_entries, self._new_lines.size+1))-1
        line_offsets = np.insert(np.cumsum(n_lines_per_entry),0, 0)
        headers = data[new_entries, 1:]
        mask  = np.ones(len(data), dtype=bool)
        mask[new_entries] = False
        sequence_lines = data[mask]
        seq_lens = sequence_lines._shape.ends[line_offsets[1:]-1]-sequence_lines._shape.starts[line_offsets[:-1]]
        sequences = RaggedArray(sequence_lines.ravel(), seq_lens)
        return SequenceEntry(headers, sequences)

    def _validate(self):
        self._is_validated = True

    @classmethod
    def from_data(cls, entries):
        name_lengths = entries.name.lengths
        sequence_lengths = entries.sequence.lengths
        n_lines = (sequence_lengths-1) // (cls.n_characters_per_line) + 1
        last_length = (sequence_lengths-1) % cls.n_characters_per_line + 1
        line_lengths = np.full(np.sum(n_lines) + n_lines.size, cls.n_characters_per_line + 1, dtype=int)
        entry_starts = np.insert(np.cumsum(n_lines+1), 0, 0)
        line_lengths[entry_starts[:-1]] = name_lengths + 2
        line_lengths[entry_starts[1:]-1] = last_length + 1
        lines = EncodedRaggedArray(
            EncodedArray(np.zeros(line_lengths.sum(), dtype=np.uint8), BaseEncoding), line_lengths)
        lines[entry_starts[:-1],1:-1] = encoded_array_from_nparray(entries.name)
        lines[entry_starts[:-1], 0] = cls._new_entry_marker
        idxs = np.delete(np.arange(len(lines)), entry_starts[:-1])
        decoded = EncodedArray(entries.sequence.encoding.decode(entries.sequence.ravel()),
                               BaseEncoding)
        lines[idxs,:-1] = EncodedRaggedArray(decoded, line_lengths[idxs]-1)
        lines[:, -1] = "\n"
        return lines.ravel()

    @classmethod
    def from_raw_buffer(cls, chunk, header_data=None):
        assert header_data is None, header_data
        chunk = EncodedArray(chunk, BaseEncoding)
        assert chunk[0] == cls._new_entry_marker, str(chunk[:100])
        new_lines = np.flatnonzero(chunk[:-1] == "\n")
        new_entries = np.flatnonzero(chunk[new_lines+1] == cls._new_entry_marker)
        if new_entries.size == 0:
            raise RuntimeError(f"No complete entry found in {cls.__name__}: \n{chunk.to_string()}\n. This can be due to badly formatted file, or because the buffer_size ({chunk.size}) is too low. Try increasing buffer_size")
        entry_starts = new_lines[new_entries]+1
        cut_chunk = chunk[:entry_starts[-1]]
        return cls(cut_chunk,
                   new_lines[:new_entries[-1]],
                   new_entries[:-1])

    def _modify_ends_for_carriage_returns(self, line_ends, data):
        if np.any((data[line_ends[:10]-1]) == '\r'):
            return line_ends-(data[line_ends-1] == '\r')
        return line_ends

    def count_entries(self):
        return len(self._new_entries)


@bnpdataclass
class FastaIdx:
    chromosome: str
    length: int
    start: int
    characters_per_line: int
    line_length: int


@bnpdataclass
class FastaIdxBuilder(FastaIdx):
    byte_size: int


class FastaIdxBuffer(MultiLineFastaBuffer):
    dataclass = FastaIdxBuilder
    
    def get_data(self):
        self.validate_if_not()
        line_starts = np.insert(self._new_lines + 1, 0, 0)
        line_ends = np.append(self._new_lines, self._data.size-1)
        entry_ends = line_ends
        line_ends = self._modify_ends_for_carriage_returns(line_ends, self._data)
        data = self._move_intervals_to_ragged_array(line_starts, line_ends)
        new_entries = np.insert(self._new_entries+1, 0, 0)
        n_lines_per_entry = np.diff(np.append(new_entries, self._new_lines.size+1))-1
        line_offsets = np.insert(np.cumsum(n_lines_per_entry), 0, 0)
        headers = data[new_entries, 1:]
        mask = np.ones(len(data), dtype=bool)
        mask[new_entries] = False
        sequence_lines = data[mask]
        ends = np.cumsum(sequence_lines.shape[-1])
        starts = np.insert(ends, 0, 0)[:-1]
        seq_lens = ends[line_offsets[1:]-1]-starts[line_offsets[:-1]]
        sequences = RaggedArray(sequence_lines.ravel(), seq_lens)

        seq_starts = line_starts[new_entries+1]
        seq_line_ends = line_ends[new_entries+1]
        chars_per_line = seq_line_ends-seq_starts
        line_lens = entry_ends[new_entries+1]-seq_starts
        return FastaIdxBuilder(headers,
                               sequences.lengths,
                               seq_starts,
                               chars_per_line,
                               line_lens+1,
                               [self._data.size]*len(headers))
