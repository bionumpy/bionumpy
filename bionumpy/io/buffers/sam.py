import numpy as np
from npstructures import RaggedArray

from ..file_buffers import TextThroughputExtractor
from ...datatypes import SAMEntry
from ..delimited_buffers import DelimitedBuffer

class SAMBufferExctractor(TextThroughputExtractor):
    def get_field_by_number(self, field_nr: int, keep_sep=False):
        if field_nr == 11:
            return self._get_extra_field()
        return super().get_field_by_number(field_nr, keep_sep)

    def _get_extra_field(self):
        starts = self._field_starts[:, -1]+self._field_lens[:, -1]+1
        lens = np.maximum(self._entry_ends-starts-1, 0)
        return self._extract_data(lens, starts)

class SAMBuffer(DelimitedBuffer):
    dataclass = SAMEntry
    COMMENT = "@"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._n_fields = 12

    @classmethod
    def _get_buffer_extractor(cls, data, delimiters, n_fields) -> TextThroughputExtractor:
        common_fields = 11# n_fields.min()
        starts = RaggedArray(delimiters[:-1] + 1, n_fields)[:, :common_fields]
        all_ends = RaggedArray(delimiters[1:], n_fields)
        all_ends = cls._modify_for_carriage_return(all_ends, data)
        ends = all_ends[:, :common_fields]
        starts, ends = (a.ravel().reshape(-1, common_fields) for a in (starts, ends))
        entry_starts = starts[:, 0]
        entry_ends = all_ends[:, -1] + 1
        return SAMBufferExctractor(data, starts, field_ends=ends, entry_starts=entry_starts, entry_ends=entry_ends)

    @classmethod
    def _get_n_fields(cls, entry_ends):
        return np.insert(np.diff(entry_ends), 0, entry_ends[0]+1)