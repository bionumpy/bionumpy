from ..encoded_array import EncodedArray, EncodedRaggedArray
from .kmers import KmerEncoder
from .rollable import RollableFunction
from ..encodings import AlphabetEncoding
from ..util import is_subclass_or_instance


class Minimizers(RollableFunction):
    def __init__(self, n_kmers, kmer_encoding=KmerEncoder):#, encoding=ACTGEncoding):
        self._n_kmers = n_kmers
        self._kmer_encoding = kmer_encoding
        self.window_size = n_kmers + kmer_encoding.window_size - 1
        self._encoding = kmer_encoding._encoding

    def __call__(self, sequence):
        kmer_hashes = self._kmer_encoding.rolling_window(sequence)
        return EncodedArray(kmer_hashes.raw().min(axis=-1), kmer_hashes.encoding)


def get_minimizers(sequence: EncodedRaggedArray, k: int, window_size: int) -> EncodedRaggedArray:
    """
    Get minimizers for sequences.
    Sequences should be encoded with an AlphabetEncoding (e.g. DNAEncoding).

    Parameters
    ----------
    sequence : EncodedRaggedArray
        Sequences to get minimizers from
    k : int
        The kmer size
    window_size : int
        The window size

    Returns
    -------
    EncodedRaggedArray
        Minimizers from the sequences.

    Examples
    --------
    >>> import bionumpy as bnp
    >>> sequences = bnp.encoded_array.as_encoded_array(["ACTG", "AAA", "TTGGC"], bnp.DNAEncoding)
    >>> bnp.sequence.get_minimizers(sequences, 2, 4)
    encoded_ragged_array([[AC],
                          [],
                          [GG, GC]], 2merEncoding(AlphabetEncoding('ACGT')))
    """
    assert is_subclass_or_instance(sequence.encoding, AlphabetEncoding), \
        "Sequence needs to be encoded with an AlphabetEncoding, e.g. DNAEncoding"
    assert k <= window_size, "kmer size must be smaller than window size"

    result = Minimizers(window_size-k+1, KmerEncoder(k, sequence.encoding)).rolling_window(sequence)
    #KmerEncoder(k, sequence.encoding).rolling_window(sequence)
    return result
