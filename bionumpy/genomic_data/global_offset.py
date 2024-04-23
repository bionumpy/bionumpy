import dataclasses
from typing import Tuple

import numpy as np
from ..encodings.string_encodings import StringEncoding
from ..encoded_array import EncodedArray, as_encoded_array
from ..bnpdataclass import replace

global_encoding = StringEncoding(["global"])


class GlobalOffset:
    def __init__(self, sequence_sizes, string_encoding=None):
        if isinstance(sequence_sizes, dict):
            self._names = as_encoded_array(list(sequence_sizes.keys()))
            self._sizes = np.array(list(sequence_sizes.values()), dtype=int)
        else:
            self._names = sequence_sizes.name
            self._sizes = sequence_sizes.size

        self._offset = np.insert(np.cumsum(self._sizes), 0, 0)
        self._old_encoding = string_encoding
        if string_encoding is None:
            self._old_encoding = StringEncoding(self._names)

    def total_size(self):
        return self._sizes.sum()

    def names(self):
        return self._names.tolist()

    def get_offset(self, seq_name):
        seq_name = as_encoded_array(seq_name, target_encoding=self._old_encoding)
        return self._offset[seq_name.raw()]

    def get_size(self, seq_name):
        seq_name = as_encoded_array(seq_name, target_encoding=self._old_encoding)
        return self._sizes[seq_name.raw()]

    def from_local_coordinates(self, sequence_name, local_offset):
        mask = local_offset >= self.get_size(sequence_name)
        if np.any(np.atleast_1d(mask)):
            raise Exception('Coordinate outside of reference:', local_offset[mask], self.get_size(sequence_name)[mask])
        return self.get_offset(sequence_name) + local_offset

    def to_local_interval(self, global_interval):
        chromosome_idxs = np.searchsorted(self._offset, global_interval.start, side="right") - 1
        start = global_interval.start - self._offset[chromosome_idxs]
        stop = global_interval.stop - self._offset[chromosome_idxs]
        assert np.all(stop <= self._sizes[chromosome_idxs])
        chromosome = EncodedArray(chromosome_idxs, self._old_encoding)
        return replace(global_interval,
                       chromosome=chromosome,
                       start=start,
                       stop=stop)

    def to_local_coordinates(self, global_offset) -> Tuple[EncodedArray, np.ndarray]:
        chromosome_idxs = np.searchsorted(self._offset, global_offset, side="right") - 1
        local_offset = global_offset - self._offset[chromosome_idxs]
        return EncodedArray(chromosome_idxs, self._old_encoding), local_offset

    def from_local_interval(self, interval, do_clip=False):
        start_offsets, stop_offsets = self.start_ends_from_intervals(interval, do_clip)
        zeros = EncodedArray(np.broadcast_to(np.array(0, dtype=np.uint8), (len(interval),)), global_encoding)
        return replace(
            interval,
            chromosome=zeros,
            start=start_offsets,
            stop=stop_offsets)

    def start_ends_from_intervals(self, interval, do_clip=False):
        chromosome = as_encoded_array(interval.chromosome, target_encoding=self._old_encoding)
        offsets = self.get_offset(chromosome)
        sizes = self.get_size(chromosome)
        if np.any(interval.start >= sizes):
            mask = interval.start >= sizes
            raise Exception(f'Starts larger than size {interval.start[mask]}, {sizes[mask]}')
        stop = interval.stop
        if do_clip:
            stop = np.minimum(stop, sizes)
        else:
            assert np.all(stop <= sizes)
        start_offsets = interval.start + offsets
        stop_offsets = stop + offsets
        return start_offsets, stop_offsets
