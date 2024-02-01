from typing import List

import numpy as np
from npstructures import RaggedArray
from npstructures.raggedshape import RaggedView2

from ..encoded_array import EncodedRaggedArray, as_encoded_array
from .exceptions import FormatException
from .file_buffers import TextBufferExtractor


class NamedBufferExtractor(TextBufferExtractor):
    def __init__(self, data, field_starts, field_lens, names):
        super().__init__(data, field_starts, field_lens=field_lens)
        self._names = names

    @classmethod
    def concatenate(cls, buffers: List['NamedBufferExtractor']):
        sizes = np.array([b._data.size for b in buffers])
        offsets = np.insert(np.cumsum(sizes), 0, 0)
        data = np.concatenate([b._data for b in buffers])
        starts = np.concatenate([b._field_starts + offset for b, offset in zip(buffers, offsets)])
        lens = np.concatenate([b._field_lens for b in buffers])
        return cls(data, starts, field_lens=lens, names=buffers[0]._names)

    def __getitem__(self, idx):
        return self.__class__(self._data, field_starts=self._field_starts[idx], field_lens=self._field_lens[idx], names=self._names)

    def get_field_by_number(self, number, keep_sep=False):
        name = self._names[number]
        return self.get_field_by_name(name, keep_sep=keep_sep)

    def has_field_number(self, number):
        name = self._names[number]
        return self.has_field_name(name)

    def has_field_name(self, name):
        starts = self._field_starts.ravel()
        mask = self._field_lens.ravel() == len(name)
        if np.any(mask):
            array = RaggedArray(self._data,
                                RaggedView2(starts[mask], np.full(mask.sum(), len(name)))).to_numpy_array()
            mask[mask] = (array == name).all(axis=-1)
        return RaggedArray(mask, self._field_starts.shape).any(axis=1)

    def get_field_by_name(self, name, keep_sep=False):
        assert name in self._names, (name, self._names)
        mask = self.has_field_mask(name)
        n_entries = len(self._field_starts)
        if not np.any(mask):
            if keep_sep:
                return EncodedRaggedArray(as_encoded_array(';'*n_entries), np.ones(n_entries, dtype=int))
            return EncodedRaggedArray(as_encoded_array(''), np.zeros(n_entries, dtype=int))
        reshaped_mask = RaggedArray(mask, self._field_starts.shape)
        line_sums = reshaped_mask.sum(axis=-1)
        if np.any(line_sums > 1):
            raise FormatException(f"Field: {name} found multiple times in buffer", line_number=np.flatnonzero(line_sums > 1)[0])
        present_mask = reshaped_mask.any(axis=-1)#line_sums.astype(bool)#

        field_starts = self._field_starts.ravel()[mask] + len(name) + 1
        lens = self._field_lens.ravel()[mask]-len(name)-1# - field_starts
        if keep_sep:
            lens += 1
        starts = np.zeros(n_entries, dtype=int)
        starts[present_mask] = field_starts
        starts = np.maximum.accumulate(starts)
        #starts = np.maximum.accumulate(np.where(present_mask, field_starts, 0))
        all_lens = np.zeros(n_entries, dtype=int)
        all_lens[present_mask] = lens
        text = EncodedRaggedArray(self._data, RaggedView2(starts, all_lens))
        return text

    def has_field_mask(self, name):
        line_len = len(name) + 1
        starts = self._field_starts.ravel()
        n_ignored_fields = 0
        while starts[len(starts) - n_ignored_fields - 1] + line_len >= self._data.size:
            n_ignored_fields += 1
        starts = starts[:len(starts) - n_ignored_fields]
        e = EncodedRaggedArray(
            self._data,
            RaggedView2(starts, [line_len] * len(starts)))
        flat_e = e.ravel()
        mask = (flat_e.reshape(-1, line_len) == name + '=').all(axis=-1)
        if n_ignored_fields:
            mask = np.append(mask, np.full(n_ignored_fields, False))
        return mask
