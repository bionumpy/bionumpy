import numpy as np
from ..io.motifs import Motif
from .rollable import RollableFunction
from ..encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from ..encodings import AlphabetEncoding
from npstructures import RaggedArray


class PositionWeightMatrix(RollableFunction):
    def __init__(self, pwm):
        self._pwm = pwm
        self.window_size = pwm.window_size

    def __call__(self, sequence: EncodedArray) -> float:
        return self._pwm.calculate_score(sequence)

        if self._encoding is not None:
            sequence = as_encoded_array(sequence, self._encoding).raw()
        scores = self._matrix[sequence, self._indices]
        return scores.sum(axis=-1)


def _pwm_from_counts(count_matrix):
    with_pseudo = count_matrix+1
    return np.log(with_pseudo/with_pseudo.sum(axis=0, keepdims=True))


class PWM:
    def __init__(self, matrix, alphabet):
        self._matrix = matrix
        self._alphabet = alphabet
        self._encoding = AlphabetEncoding(alphabet)
        self._indices = np.arange(self.window_size)

    @property
    def alphabet(self):
        return self._alphabet

    @property
    def window_size(self):
        return self._matrix.shape[-1]

    def calculate_score(self, sequence):
        sequence = as_encoded_array(sequence, self._encoding)
        assert sequence.encoding == self._encoding
        assert sequence.shape[-1] == self.window_size
        scores = self._matrix[sequence.raw(), self._indices]
        return scores.sum(axis=-1)

    @classmethod
    def from_dict(cls, dictionary: dict):
        alphabet = "".join(dictionary.keys())
        matrix = np.array(list(dictionary.values()))
        return cls(matrix, alphabet)

    @classmethod
    def from_motif(cls, motif: Motif):
        return cls(motif.matrix, motif.alphabet)

    @classmethod
    def from_counts(cls, motif: Motif):
        return cls(_pwm_from_counts(motif.matrix), motif.alphabet)



def get_motif_scores(sequence: EncodedRaggedArray, pwm: PWM) -> RaggedArray:
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
    pwm = PositionWeightMatrix(pwm)
    return pwm.rolling_window(sequence)
