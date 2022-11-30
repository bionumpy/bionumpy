import numpy as np
from bionumpy.rollable import RollableFunction
from bionumpy.encoded_array import EncodedArray, EncodedRaggedArray
from bionumpy import as_encoded_array
from bionumpy.encodings import DNAEncoding, AlphabetEncoding
from npstructures import RaggedArray


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


def _pwm_from_counts(count_matrix):
    with_pseudo = count_matrix+1
    return np.log(with_pseudo/with_pseudo.sum(axis=0, keepdims=True))


def get_motif_scores(sequence: EncodedRaggedArray, motif: PositionWeightMatrix) -> RaggedArray:
    """Computes motif scores for a motif on a sequence.
    Returns a RaggedArray with the score at each position in every read.

    Parameters
    ----------
    sequence: EncodedRaggedArray
    motif: PositionWeightMatrix

    Returns
    -------
    RaggedArray
        A numeric RaggedArray. Contains one row for every read
        with the scores for every position of that read.

    Examples
    --------

    """
    pwm = PositionWeightMatrix(_pwm_from_counts(motif.matrix), AlphabetEncoding(motif.alphabet))
    return pwm.rolling_window(sequence)
