from ..file_buffers import TextThroughputExtractor
from ...datatypes import SAMEntry
from ..delimited_buffers import DelimitedBuffer


class SAMBuffer(DelimitedBuffer):
    dataclass = SAMEntry
    COMMENT = "@"

    @classmethod
    def _get_buffer_extractor(cls, data, delimiters, n_cols) -> TextThroughputExtractor:
        starts = delimiters[:-1].reshape(-1, n_cols) + 1
        ends = delimiters[1:].reshape(-1, n_cols)
        ends = cls._modify_for_carriage_return(ends, data)
        entry_starts = starts[:, 0]
        entry_ends = ends[:, -1] + 1
        return TextThroughputExtractor(data, starts, field_ends=ends, entry_starts=entry_starts, entry_ends=entry_ends)


