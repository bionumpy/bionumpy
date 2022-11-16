import numpy as np
from bionumpy.rollable import RollableFunction
from bionumpy.encoded_array import as_encoded_array, EncodedArray
from bionumpy.encodings import DNAEncoding, AlphabetEncoding


class PositionWeightMatrix(RollableFunction):
    def __init__(self, matrix, encoding=None):
        self._matrix = np.asanyarray(matrix)
        self.window_size = self._matrix.shape[-1]
        self._indices = np.arange(self.window_size)
        self._encoding = encoding

    def __call__(self, sequence: EncodedArray) -> float:
        if self._encoding is not None:
            sequence = as_encoded_array(sequence, self._encoding).raw()
        scores = self._matrix[sequence, self._indices]
        return scores.sum(axis=-1)


def pwm_from_counts(count_matrix):
    with_pseudo = count_matrix+1
    return np.log(with_pseudo/with_pseudo.sum(axis=0, keepdims=True))


def get_motif_scores(sequence, motif):
    pwm = PositionWeightMatrix(pwm_from_counts(motif.matrix), AlphabetEncoding(motif.alphabet))
    return pwm.rolling_window(sequence)
