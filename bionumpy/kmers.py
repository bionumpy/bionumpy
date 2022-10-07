import numpy as np
from .rollable import RollableFunction
from .encodings import ACTGEncoding
from .sequences import Sequence, as_encoded_sequence_array
from npstructures.bitarray import BitArray
from .util import convolution
import logging
logger = logging.getLogger(__name__)


class KmerEncoding(RollableFunction):

    def __init__(self, k, alphabet_encoding):
        self.window_size = k
        self._k = k
        self._encoding = alphabet_encoding
        self._alphabet_size = alphabet_encoding.encoding.alphabet_size
        self._convolution = self._alphabet_size ** np.arange(self._k)

    def __call__(self, sequence: Sequence) -> int:
        sequence = as_encoded_sequence_array(sequence, encoding=self._encoding)
        return np.asarray(sequence).dot(self._convolution)

    def inverse(self, kmer_hash: int) -> Sequence:
        return ((kmer_hash[:, np.newaxis] // self._convolution) % self._alphabet_size).view(self._encoding)
        # s.encoding=self._encoding

    def sample_domain(self, n):
        return (np.random.randint(0, self._alphabet_size, size=self._k * n).reshape(
            n, self._k)).view(self._encoding)


@convolution
def fast_hash(sequence, k, encoding=None):
    sequence = as_encoded_sequence_array(sequence, ACTGEncoding)
    if encoding:
        sequence = encoding.encode(sequence)

    bit_array = BitArray.pack(sequence, bit_stride=2)
    hashes = bit_array.sliding_window(k)
    return hashes


