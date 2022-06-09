from bionumpy.rollable import RollableFunction
from bionumpy.sequences import as_sequence_array
import numpy as np
from bionumpy.sequences import Sequence

class StringMatcher(RollableFunction):
    def __init__(self, matching_sequence, encoding):
        self._encoding = encoding
        self._matching_sequence_array = as_sequence_array(matching_sequence, encoding=encoding)

    @property
    def window_size(self):
        return len(self._matching_sequence_array)

    def __call__(self, sequence):
        return np.all(sequence == self._matching_sequence_array, axis=-1)


# class MaskedStringMatcher(RollableFunction):
#     def __init__(self, matching_sequence_array, mask, encoding):
#         assert isinstance(matching_sequence_array, Sequence)
#         self._matching_sequence_array = matching_sequence_array
#         self._mask = mask
#         self._encoding = encoding
#
#     @property
#     def window_size(self):
#         return len(self._matching_sequence_array)
#
#     def __call__(self, sequence):
#         return np.all(sequence == self._matching_sequence_array, axis=-1)
#
#
