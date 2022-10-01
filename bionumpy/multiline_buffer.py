from npstructures import RaggedArray
from .file_buffers import FileBuffer, NEWLINE
from .sequences import EncodedArray, ASCIIText
from .datatypes import SequenceEntry
import numpy as np


class MultiLineBuffer(FileBuffer):
    pass


class MultiLineFastaBuffer(MultiLineBuffer):
    _new_entry_marker = ">"
    n_characters_per_line = 80
    dataclass = SequenceEntry

    def __init__(self, data, new_lines, new_entries):
        super().__init__(data, new_lines)
        self._new_entries = new_entries
    
    def get_data(self):
        self.validate_if_not()
        line_starts = np.insert(self._new_lines + 1, 0, 0)
        line_ends = np.append(self._new_lines, self._data.size-1)
        data = self._move_intervals_to_ragged_array(line_starts, line_ends)
        new_entries = np.insert(self._new_entries+1, 0, 0)
        n_lines_per_entry = np.diff(np.append(new_entries, self._new_lines.size+1))-1
        line_offsets = np.insert(np.cumsum(n_lines_per_entry),0, 0)
        headers = data[new_entries, 1:]
        mask  = np.ones(len(data), dtype=bool)
        mask[new_entries] = False
        sequence_lines = data[mask]
        seq_lens = sequence_lines.shape.ends[line_offsets[1:]-1]-sequence_lines.shape.starts[line_offsets[:-1]]
        sequences = RaggedArray(sequence_lines.ravel(), seq_lens)
        return SequenceEntry(headers, sequences)

    def _validate(self):
        self._is_validated=True

    @classmethod
    def from_data(cls, entries):
        name_lengths = entries.name.shape.lengths
        sequence_lengths = entries.sequence.shape.lengths
        n_lines = (sequence_lengths-1) // (cls.n_characters_per_line) + 1
        last_length = (sequence_lengths-1) % cls.n_characters_per_line + 1
        line_lengths = np.full(np.sum(n_lines) + n_lines.size, cls.n_characters_per_line + 1, dtype=int)
        entry_starts = np.insert(np.cumsum(n_lines+1), 0, 0)
        line_lengths[entry_starts[:-1]] = name_lengths + 2
        line_lengths[entry_starts[1:]-1] = last_length + 1
        lines = RaggedArray(np.zeros(line_lengths.sum(), dtype=np.uint8).view(EncodedArray), line_lengths)
        lines[entry_starts[:-1],1:-1] = entries.name
        lines[entry_starts[:-1], 0] = cls._new_entry_marker
        idxs = np.delete(np.arange(len(lines)), entry_starts[:-1])
        lines[idxs,:-1] = RaggedArray(entries.sequence.ravel(), line_lengths[idxs]-1)
        lines[:, -1] = "\n"
        return lines.ravel()
        
    @classmethod
    def from_raw_buffer(cls, chunk, header_data=None):
        assert header_data is None, header_data
        chunk = chunk.view(ASCIIText)
        assert chunk[0] == cls._new_entry_marker, str(chunk[:100])
        new_lines = np.flatnonzero(chunk[:-1] == "\n")
        new_entries = np.flatnonzero(chunk[new_lines+1] == cls._new_entry_marker)
        entry_starts = new_lines[new_entries]+1
        cut_chunk = chunk[:entry_starts[-1]]
        return cls(cut_chunk,
                   new_lines[:new_entries[-1]],
                   new_entries[:-1])
