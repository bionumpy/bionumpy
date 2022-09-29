import numpy as np
from .rollable import RollableFunction
from .sequences import as_encoded_sequence_array, Sequence
from .encodings import ACTGEncoding


class PositionWeightMatrix(RollableFunction):
    def __init__(self, matrix, encoding=None):
        self._matrix = np.asanyarray(matrix)
        self.window_size = self._matrix.shape[-1]
        self._indices = np.arange(self.window_size)
        self._encoding = encoding

    def __call__(self, sequence: Sequence) -> float:
        if self._encoding is not None:
            sequence = as_encoded_sequence_array(sequence, self._encoding)
        scores = self._matrix[np.asarray(sequence), self._indices]
        return scores.sum(axis=-1)


def pwm_from_counts(count_matrix):
    with_pseudo = count_matrix+1
    return np.log(with_pseudo/with_pseudo.sum(axis=0, keepdims=True))
