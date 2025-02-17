import numpy as np
from ..streams import streamable
from .count_encoded import count_encoded, EncodedCounts
from ..encodings.exceptions import EncodingError
from ..encodings.kmer_encodings import KmerEncoding
from .rollable import RollableFunction
from bionumpy.encodings import DNAEncoding
from ..encodings.alphabet_encoding import AlphabetEncoding
from ..encoded_array import EncodedArray, EncodedRaggedArray, BaseEncoding, change_encoding
from ..encoded_array import as_encoded_array
from npstructures.bitarray import BitArray
from ..util import as_strided, is_subclass_or_instance
import logging
logger = logging.getLogger(__name__)


class KmerEncoder(RollableFunction):
    def __init__(self, k, alphabet_encoding):
        self.window_size = k
        self._k = k
        self._encoding = alphabet_encoding
        self._alphabet_size = alphabet_encoding.alphabet_size
        self._convolution = self._alphabet_size ** np.arange(self._k)

    def __call__(self, sequence: EncodedArray) -> int:
        sequence = as_encoded_array(sequence, target_encoding=self._encoding)
        return EncodedArray(sequence.data.dot(self._convolution), KmerEncoding(self._encoding, self._k))

    def inverse(self, kmer_hash: int) -> EncodedArray:
        return EncodedArray(((kmer_hash[:, np.newaxis] // self._convolution) % self._alphabet_size), self._encoding)

    def sample_domain(self, n):
        return EncodedArray((np.random.randint(0, self._alphabet_size, size=self._k * n).reshape(n, self._k)), self._encoding)


def get_kmers(sequence: EncodedRaggedArray, k: int) -> EncodedArray:
    """
    Get kmers for sequences.
    Sequences should be encoded with an AlphabetEncoding (e.g. DNAEncoding).
    Use bnp.change_encoding if your sequences do not have a suitable encoding.

    Parameters
    ----------
    sequence : EncodedRaggedArray
        Sequences to get kmers from
    k : int
        The kmer size (1-31)

    Returns
    -------
    EncodedRaggedArray
        Kmers from the sequences.

    Examples
    --------
    >>> import bionumpy as bnp
    >>> sequences = bnp.encoded_array.as_encoded_array(["ACTG", "AAA", "TTGGC"], bnp.DNAEncoding)
    >>> bnp.sequence.get_kmers(sequences, 3)
    encoded_ragged_array([[ACT, CTG],
                          [AAA],
                          [TTG, TGG, GGC]], 3merEncoding(AlphabetEncoding('ACGT')))

    >>> sequences = bnp.open("example_data/big.fq.gz").read().sequence
    >>> sequences = bnp.change_encoding(sequences, bnp.DNAEncoding)
    >>> bnp.sequence.get_kmers(sequences, 31)[0, 0:3]  # first three kmers of first sequence
    encoded_array([CGGTAGCCAGCTGCGTTCAGTATGGAAGATT, GGTAGCCAGCTGCGTTCAGTATGGAAGATTT, GTAGCCAGCTGCGTTCAGTATGGAAGATTTG], 31merEncoding(AlphabetEncoding('ACGT')))
    """

    assert 0 < k < 32, "k must be larger than 0 and smaller than 32"
    if sequence.encoding == BaseEncoding:
        try:
            sequence = change_encoding(sequence, DNAEncoding)
        except EncodingError:
            logging.error("Tried to change encoding of sequences to DNAEncoding, but failed. "
                          "Make sure your sequences are valid DNA, only containing A, C, G, and T")
            raise

    assert is_subclass_or_instance(sequence.encoding, AlphabetEncoding), \
        "Sequence needs to be encoded with an AlphabetEncoding, e.g. DNAEncoding. " \
        "Change encoding of your sequences by using e.g. bnp.change_encoding(sequences, bnp.DNAEncoding)"

    if sequence.encoding.alphabet_size == 4:
        # use the faster _get_dna_kmers
        result = _get_dna_kmers(sequence, k)
        return result

    return KmerEncoder(k, sequence.encoding).rolling_window(sequence)


def convolution(func):
    def new_func(_sequence, window_size, *args, **kwargs):
        shape, sequence = (_sequence.shape, _sequence.ravel())
        convoluted = func(sequence, window_size, *args, **kwargs)
        if not isinstance(shape[-1], np.ndarray):
            out = as_strided(convoluted, shape)
        else:
            out = EncodedRaggedArray(convoluted, shape, safe_mode=False) # FIXME: Should not use unsafe here
        if window_size == 1:
            return out
        return out[..., : (-window_size + 1)]

    return new_func


@convolution
def _get_dna_kmers(sequence, k):
    """
    A faster version alternative to get_kmers that can be used with sequences encoded
     with an AlphabetEncoding with alphabet size **exactly equal to 4**.

    This function is called by get_kmers when alphabet size is 4 to get some speedup.

    See get_kmers for documentation on use.
    """
    assert 0 < k < 32, "k must be larger than 0 and smaller than 32"
    assert isinstance(sequence, EncodedArray), sequence
    assert is_subclass_or_instance(sequence.encoding, AlphabetEncoding)
    assert sequence.encoding.alphabet_size == 4, \
        "Only supported for sequence encoded with an AlphabetEncoding with alphabetsize 4"

    bit_array = BitArray.pack(sequence.data, bit_stride=2)
    hashes = bit_array.sliding_window(k)
    assert hashes.dtype == np.uint64
    hashes = hashes.view(np.int64)
    output_encoding = KmerEncoding(sequence.encoding, k)
    return EncodedArray(hashes, output_encoding)


@streamable(sum)
def count_kmers(sequence: EncodedRaggedArray, k: int, axis=None) -> EncodedCounts:
    """Count the occurance of each kmer in a set of sequences

    For each possible kmer, count the number of occurances in a (set
    of) sequence(s)

    Parameters
    ----------
    sequence : EncodedRaggedArray
    k : int
        Length of the kmer
    axis : int

    """
    kmers = get_kmers(sequence, k)
    return count_encoded(kmers, axis=axis)
