from .file_buffers import FileBuffer, NEWLINE
from .datatypes import SequenceEntry
from .sequences import Sequences
import numpy as np

class MultiLineBuffer(FileBuffer):
    pass


class MultiLineFastaBuffer(MultiLineBuffer):
    _new_entry_marker = ord(">")

    def __init__(self, data, new_lines, new_entries):
        super().__init__(data, new_lines)
        self._new_entries = new_entries
    
    def get_data(self):
        line_starts = np.insert(self._new_lines + 1, 0, 0)
        line_ends = np.append(self._new_lines, self._data.size)
        data = self._move_intervals_to_ragged_array(line_starts, line_ends)
        new_entries = np.insert(self._new_entries+1, 0, 0)
        n_lines_per_entry = np.diff(np.append(new_entries, self._new_lines.size+1))-1
        line_offsets = np.insert(np.cumsum(n_lines_per_entry),0, 0)
        headers = data[new_entries, 1:]
        mask  = np.ones(len(data), dtype=bool)
        mask[new_entries] = False
        print(headers)
        print(n_lines_per_entry)
        sequence_lines = data[mask]
        seq_lens = sequence_lines.shape.ends[line_offsets[1:]-1]-sequence_lines.shape.starts[line_offsets[:-1]]
        sequences = Sequences(sequence_lines.ravel(), seq_lens)
        print(sequences)
        return SequenceEntry(headers, sequences)

    @classmethod
    def from_raw_buffer(cls, chunk, header_data=None):
        assert header_data is None
        assert chunk[0] == cls._new_entry_marker
        new_lines = np.flatnonzero(chunk == NEWLINE)
        new_entries = np.flatnonzero(chunk[new_lines+1] == cls._new_entry_marker)
        entry_starts = new_lines[new_entries]+1
        print(new_lines, new_entries)
        return cls(chunk[:entry_starts[-1]-1],
                   new_lines[:new_entries[-1]],
                   new_entries[:-1])
