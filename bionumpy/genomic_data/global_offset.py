import dataclasses
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
        return self.get_offset(sequence_name) + local_offset

    def to_local_interval(self, global_interval):
        chromosome_idxs = np.searchsorted(self._offset, global_interval.start, side="right")-1
        start = global_interval.start-self._offset[chromosome_idxs]
        stop = global_interval.stop-self._offset[chromosome_idxs]
        assert np.all(stop <= self._sizes[chromosome_idxs])
        chromosome = EncodedArray(chromosome_idxs, self._old_encoding)
        return dataclasses.replace(global_interval,
                                   chromosome=chromosome,
                                   start=start,
                                   stop=stop)

    def from_local_interval(self, interval, do_clip=False):
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
        return replace(
            interval,
            chromosome=EncodedArray(
                np.broadcast_to(np.array(0,dtype=np.uint8), (len(interval), )), global_encoding),
            start=interval.start+offsets,
            stop=stop+offsets)
