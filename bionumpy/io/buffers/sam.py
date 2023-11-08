import numpy as np
from npstructures import RaggedArray

from ..file_buffers import TextThroughputExtractor
from ...datatypes import SAMEntry
from ..delimited_buffers import DelimitedBuffer


class SAMBuffer(DelimitedBuffer):
    dataclass = SAMEntry
    COMMENT = "@"

    @classmethod
    def _get_buffer_extractor(cls, data, delimiters, n_fields) -> TextThroughputExtractor:
        common_fields = 11# n_fields.min()
        starts = RaggedArray(delimiters[:-1] + 1, n_fields)[:, :common_fields]
        ends= RaggedArray(delimiters[1:], n_fields)[:, :common_fields].ravel()
        starts, ends = (a.ravel().reshape(-1, common_fields) for a in (starts, ends))
        ends = cls._modify_for_carriage_return(ends, data)
        entry_starts = starts[:, 0]
        entry_ends = ends[:, -1] + 1
        return TextThroughputExtractor(data, starts, field_ends=ends, entry_starts=entry_starts, entry_ends=entry_ends)

    @classmethod
    def _get_n_fields(cls, entry_ends):
        return np.insert(np.diff(entry_ends), 0, entry_ends[0]+1)