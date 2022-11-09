import numpy as np
from .rollable import RollableFunction
from .encodings import DNAEncoding
from .encodings.alphabet_encoding import AlphabetEncoding
from .encoded_array import EncodedArray, as_encoded_array
from npstructures.bitarray import BitArray
from .util import convolution, is_subclass_or_instance
import logging
logger = logging.getLogger(__name__)


class KmerEncoding(RollableFunction):

    def __init__(self, k, alphabet_encoding):
        self.window_size = k
        self._k = k
        self._encoding = alphabet_encoding
        self._alphabet_size = alphabet_encoding.alphabet_size
        self._convolution = self._alphabet_size ** np.arange(self._k)

    def __call__(self, sequence: EncodedArray) -> int:
        sequence = as_encoded_array(sequence, target_encoding=self._encoding)
        return sequence.data.dot(self._convolution)

    def inverse(self, kmer_hash: int) -> EncodedArray:
        return EncodedArray(((kmer_hash[:, np.newaxis] // self._convolution) % self._alphabet_size), self._encoding)
        # s.encoding=self._encoding

    def sample_domain(self, n):
        return EncodedArray((np.random.randint(0, self._alphabet_size, size=self._k * n).reshape(n, self._k)), self._encoding)


@convolution
def fast_hash(sequence, k, encoding=None):
    assert isinstance(sequence, EncodedArray), sequence
    assert is_subclass_or_instance(sequence.encoding, AlphabetEncoding)

    bit_array = BitArray.pack(sequence.data, bit_stride=2)
    hashes = bit_array.sliding_window(k)
    return hashes


