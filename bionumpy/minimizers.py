from . import EncodedArray
from .kmers import KmerEncoder, get_kmers
from .rollable import RollableFunction
from .encodings import DNAEncoding, AlphabetEncoding
from .util import is_subclass_or_instance


class Minimizers(RollableFunction):
    def __init__(self, n_kmers, kmer_encoding=KmerEncoder):#, encoding=ACTGEncoding):
        self._n_kmers = n_kmers
        self._kmer_encoding = kmer_encoding
        self.window_size = n_kmers + kmer_encoding.window_size - 1
        # self._encoding = encoding

    def __call__(self, sequence):
        kmer_hashes = self._kmer_encoding.rolling_window(sequence)
        return EncodedArray(kmer_hashes.raw().min(axis=-1), kmer_hashes.encoding)


def get_minimizers(sequence, k, window_size):
    assert is_subclass_or_instance(sequence.encoding, AlphabetEncoding), \
        "Sequence needs to be encoded with an AlphabetEncoding, e.g. DNAEncoding"

    result = Minimizers(window_size-k+1, KmerEncoder(k, sequence.encoding)).rolling_window(sequence)
    #KmerEncoder(k, sequence.encoding).rolling_window(sequence)
    return result
