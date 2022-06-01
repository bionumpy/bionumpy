import numpy as np
from .rollable import RollableFunction
# from .encodings import ACGTEncoding
LobProbabilities = None


class PositionWeightMatrix(RollableFunction):
    def __init__(self, matrix):
        self._matrix = np.asanyarray(matrix)
        self.window_size = self._matrix.shape[-1]
        self._indices = np.arange(self.window_size)

    def __call__(self, sequence):
        scores = self._matrix[sequence, self._indices]
        return scores.sum(axis=-1)

    @property
    def range(self):
        return (-np.inf, 0)

    def sample_domain(self, n):
        return self.encoding.sample_domain(n*self.window_size).reshape(-1, self.window_size)
