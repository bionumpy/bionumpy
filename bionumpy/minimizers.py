from .kmers import KmerEncoding
from .rollable import RollableFunction
from .encodings import ACTGEncoding


class Minimizers(RollableFunction):
    def __init__(self, n_kmers, kmer_encoding=KmerEncoding):#, encoding=ACTGEncoding):
        self._n_kmers = n_kmers
        self._kmer_encoding = kmer_encoding
        self.window_size = n_kmers + kmer_encoding.window_size - 1
        # self._encoding = encoding

    def __call__(self, sequence):
        kmer_hashes = self._kmer_encoding.rolling_window(sequence)
        return kmer_hashes.min(axis=-1)
