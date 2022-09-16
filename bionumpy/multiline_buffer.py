from file_buffer import FileBuffer


class MultiLineBuffer(FileBuffer):
    pass


class MultiLineFastaBuffer:

    @classmethod
    def from_raw_buffer(cls, chunk, header_data=None):
        assert header_data is None
        new_lines = np.flatnonzero(chunk == NEWLINE)
        new_entries = np.flatnonzero(chunk(new_lines+1)) ==

        n_lines = new_lines.size
        assert n_lines >= cls.n_lines_per_entry, "No complete entry in buffer"
        new_lines = new_lines[: n_lines - (n_lines % cls.n_lines_per_entry)]
        return cls(chunk[: new_lines[-1] + 1], new_lines)
        pass
