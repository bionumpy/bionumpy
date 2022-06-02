import numpy as np
from .rollable import RollableFunction
from .encodings import ACTGEncoding
from .sequences import Sequence, as_sequence_array
from npstructures.bitarray import BitArray
from .util import convolution
import logging
logger = logging.getLogger(__name__)


class KmerEncoding(RollableFunction):

    def __init__(self, k, alphabet_size=4, encoding=ACTGEncoding):
        self.window_size = k
        self._k = k
        self._alphabet_size = alphabet_size
        self._convolution = self._alphabet_size ** np.arange(self._k)
        self._encoding = encoding

    def __call__(self, sequence: Sequence) -> int:
        sequence = as_sequence_array(sequence, encoding=self._encoding)
        return sequence.dot(self._convolution)

    def inverse(self, array: int) -> Sequence:
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
