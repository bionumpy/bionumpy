import numpy as np
from npstructures import RaggedArray

from abc import abstractmethod


class RollableFunction:

    @abstractmethod
    def __call__(self, sequence):
        return NotImplemented

    def rolling_window(self, _sequence, window_size=None):

        if window_size is None:
            window_size = self.window_size


        shape, sequence = (_sequence.shape, _sequence.ravel())
        windows = np.lib.stride_tricks.sliding_window_view(
            sequence, window_size)
        convoluted = self(windows)
        if isinstance(_sequence, RaggedArray):
            out = RaggedArray(convoluted, shape)
        elif isinstance(_sequence, np.ndarray):
            out = np.lib.stride_tricks.as_strided(convoluted, shape)
        return out[..., :(-window_size+1)]
