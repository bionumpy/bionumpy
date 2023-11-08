import numpy as np

from ..file_buffers import TextThroughputExtractor
from ...datatypes import SAMEntry
from ..delimited_buffers import DelimitedBuffer


class SAMBuffer(DelimitedBuffer):
    dataclass = SAMEntry
    COMMENT = "@"

    @classmethod
    def __get_buffer_extractor(cls, data, delimiters, n_cols) -> TextThroughputExtractor:
        starts = delimiters[:-1].reshape(-1, n_cols) + 1
        ends = delimiters[1:].reshape(-1, n_cols)
        ends = cls._modify_for_carriage_return(ends, data)
        entry_starts = starts[:, 0]
        entry_ends = ends[:, -1] + 1
        return TextThroughputExtractor(data, starts, field_ends=ends, entry_starts=entry_starts, entry_ends=entry_ends)

    @classmethod
    def __from_raw_buffer(cls, chunk: np.ndarray, header_data=None) -> "DelimitedBuffer":
        chunk = EncodedArray(chunk, BaseEncoding)
        mask = chunk == NEWLINE
        mask |= chunk == cls.DELIMITER
        delimiters = np.flatnonzero(mask)
        n_fields = np.flatnonzero(chunk[delimiters] == '\n')
        if n_fields.size == 0:
            logging.warning("Foud no new lines. Chunk size may be too low. Try increasing")
            raise
        n_fields = n_fields[0] + 1
        new_lines = delimiters[(n_fields - 1)::n_fields]
        delimiters = np.concatenate(([-1], delimiters[:n_fields * len(new_lines)]))
        buffer_extractor = cls._get_buffer_extractor(
            chunk[:new_lines[-1] + 1], delimiters, n_fields)
        return cls(buffer_extractor, header_data)
