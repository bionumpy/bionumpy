import dataclasses
import numpy as np
from ..encodings.string_encodings import StringEncoding
from ..encoded_array import EncodedArray, as_encoded_array

global_encoding = StringEncoding(["global"])


class GlobalOffset:
    def __init__(self, sequence_sizes):
        self._names = sequence_sizes.name
        self._sizes = sequence_sizes.size
        self._offset = np.insert(np.cumsum(self._sizes), 0, 0)
        self._old_encoding = StringEncoding(self._names)

    def get_offset(self, seq_name):
        seq_name = as_encoded_array(seq_name, target_encoding=self._old_encoding)
        return self._offset[seq_name.raw()]

    def get_size(self, seq_name):
        seq_name = as_encoded_array(seq_name, target_encoding=self._old_encoding)
        return self._sizes[seq_name.raw()]

    def from_local_coordinates(self, sequence_name, local_offset):
        return self.get_offset(sequence_name) + local_offset

    def from_local_interval(self, interval, do_clip=False):
        offsets = self.get_offset(interval.chromosome)
        sizes = self.get_size(interval.chromosome)
        assert np.all(interval.start < sizes)
        stop = interval.stop
        if do_clip:
            stop = np.minimum(stop, sizes)
        else:
            assert np.all(stop <= sizes)
        return dataclasses.replace(
            interval,
            chromosome=EncodedArray(np.zeros(len(interval), dtype=np.uint8), global_encoding),
            start=interval.start+offsets,
            stop=stop+offsets)
