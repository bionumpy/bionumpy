import numpy as np
from .rollable import RollableFunction
from .encodings import ACTGEncoding
from .sequences import Sequence
from npstructures.bitarray import BitArray
from .util import convolution
import logging
logger = logging.getLogger(__name__)


class KmerEncoding(RollableFunction):

    def __init__(self, k, alphabet_size=4):
        self.window_size = k
        self._k = k
        self._alphabet_size = alphabet_size
        self._convolution = self._alphabet_size ** np.arange(self._k)

    def __call__(self, sequence: Sequence) -> np.ndarray:
        return sequence.dot(self._convolution)

    def inverse(self, array: np.ndarray) -> Sequence:
        return Sequence.from_array(
            (array[:, np.newaxis] // self._convolution) % self._alphabet_size)

    def sample_domain(self, n):
        return np.random.randint(0, self._alphabet_size, size=self._k * n).reshape(
            n, self._k
        )


@convolution
def fast_hash(sequence, k):
    encoded = ACTGEncoding.from_bytes(sequence)
    bit_array = BitArray.pack(encoded, bit_stride=2)
    hashes = bit_array.sliding_window(k)
    return hashes
