import numpy as np
from .rollable import RollableFunction


class PositionWeightMatrix(RollableFunction):
    def __init__(self, matrix):
        self._matrix = np.asanyarray(matrix)
        self.window_size = self._matrix.shape[-1]
        self._indices = np.arange(self.window_size)

    def __call__(self, sequence):
        scores = self._matrix[sequence, self._indices]
        return scores.sum(axis=-1)
