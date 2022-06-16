import numpy as np
from .rollable import RollableFunction
from .encodings import ACTGEncoding
from .sequences import Sequence, as_sequence_array
from npstructures.bitarray import BitArray
from .util import convolution
import logging
logger = logging.getLogger(__name__)


class KmerEncoding(RollableFunction):

<<<<<<< HEAD
    def __init__(self, k, alphabet_encoding=None, alphabet_size=None):
        self.window_size = k
=======
    def hash(self, kmer):
        return kmer.dot(self._powers)


@convolution
def fast_hash(sequence, k, do_encoding=True):
    if do_encoding:
        encoded = ACTGEncoding.from_bytes(sequence)
    else:
        encoded = sequence
    bit_array = BitArray.pack(encoded, bit_stride=2)
    hashes = bit_array.sliding_window(k)
    return hashes


def hash_sequences(sequences, k, hash_func):
    sequence = sequences.ravel()
    kmers = hash_func(sequence, k)
    ragged_array = RaggedArray(kmers, sequences.shape)
    return ragged_array[:, :-k+1]


class KmerEncoding:
    def __init__(self, k, alphabet_size=4):
>>>>>>> fast-hash-encoding-fix
        self._k = k
        self._encoding = alphabet_encoding
        if alphabet_encoding is not None:
            alphabet_size = alphabet_encoding.alphabet_size
        self._alphabet_size = alphabet_size
        self._convolution = self._alphabet_size ** np.arange(self._k)

    def __call__(self, sequence: Sequence) -> int:
        if self._encoding is not None:
            sequence = as_sequence_array(sequence, encoding=self._encoding)
        return sequence.dot(self._convolution)

    def inverse(self, kmer_hash: int) -> Sequence:
        return Sequence.from_array(
            (kmer_hash[:, np.newaxis] // self._convolution) % self._alphabet_size)

    def sample_domain(self, n):
        s = Sequence.from_array(np.random.randint(0, self._alphabet_size, size=self._k * n).reshape(
            n, self._k
        ))
        s.encoding = self._encoding
        return s


@convolution
def fast_hash(sequence, k):
    encoded = ACTGEncoding.encode(sequence)
    bit_array = BitArray.pack(encoded, bit_stride=2)
    hashes = bit_array.sliding_window(k)
    return hashes
