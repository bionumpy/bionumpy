import numpy as np
import typing
from ..io.motifs import Motif
from .rollable import RollableFunction
from ..encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from ..encodings import AlphabetEncoding
from npstructures import RaggedArray


class PositionWeightMatrix(RollableFunction):
    def __init__(self, pwm):
        self._pwm = pwm
        self._encoding = pwm._encoding
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

    def calculate_scores(self, sequence: EncodedArray):
        sequence = as_encoded_array(sequence, self._encoding)
        assert sequence.encoding == self._encoding
        scores = np.zeros(sequence.size, dtype=float)
        m = self._matrix.T.copy()
        for offset, row in enumerate(m):
            scores[:scores.size-offset] += row[sequence[offset:].raw()]
        return scores

    @classmethod
    def from_dict(cls, dictionary: dict):
        alphabet = "".join(dictionary.keys())
        matrix = np.log(np.array(list(dictionary.values())))
        return cls(matrix, alphabet)

    @classmethod
    def from_motif(cls, motif: Motif):
        return cls(motif.matrix, motif.alphabet)

    @classmethod
    def from_counts(cls, counts: typing.Union[Motif, dict]):
        if isinstance(counts, Motif):
            return cls(_pwm_from_counts(counts.matrix), counts.alphabet)
        else:
            return cls(_pwm_from_counts(np.array(counts.values()), "".join(counts.keys())))


def get_motif_scores_old(sequence: EncodedRaggedArray, pwm: PWM) -> RaggedArray:
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
    sequence = as_encoded_array(sequence)
    flat_sequence, shape = (sequence.ravel(), sequence.shape)
    scores = pwm.calculate_scores(flat_sequence)
    if isinstance(sequence, EncodedRaggedArray):
        scores = RaggedArray(scores, shape[-1])
    return scores[..., :(-pwm.window_size+1)]
