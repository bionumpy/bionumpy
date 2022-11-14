from . import EncodedRaggedArray
from .kmers import KmerEncoder
import numpy as np
from collections import defaultdict


class KmerIndex:
    def __init__(self, k, lookup, sequences_encoding):
        self._k = k
        self._lookup = lookup
        self._sequences_encoding = sequences_encoding

    @classmethod
    def create_index(cls, sequences: EncodedRaggedArray, k: int) -> "KmerIndex":
        func = KmerEncoder(k=k, alphabet_encoding=sequences.encoding).rolling_window
        kmers = func(sequences).raw()
        unique_kmers = np.unique(kmers.ravel())
        lookup = defaultdict(list)
        for kmer in unique_kmers:
            lookup[int(kmer)] = cls._get_index_for_kmer(kmers, kmer)
        return cls(k, lookup, sequences.encoding)

    @classmethod
    def _get_index_for_kmer(cls, kmers, kmer):
        return np.flatnonzero(np.any(kmers == kmer, axis=-1))

    def get_kmer_indices(self, kmer):
        if isinstance(kmer, str):
            assert len(kmer) == self._k
            encoded_kmer = KmerEncoder(self._k, self._sequences_encoding)(kmer).raw()
            return self._lookup[int(encoded_kmer)]
        else:
            return self._lookup[int(kmer)]


class KmerLookup:
    def __init__(self, kmer_index: KmerIndex, sequences: EncodedRaggedArray):
        self._kmer_index = kmer_index
        self._sequences = sequences

    @classmethod
    def create_lookup(cls, sequences: EncodedRaggedArray, k: int) -> "KmerLookup":
        index = KmerIndex.create_index(sequences=sequences, k=k)
        return cls(index, sequences)

    def get_sequences_with_kmer(self, kmer):
        kmer_indices = self._kmer_index.get_kmer_indices(kmer)
        return self._sequences[kmer_indices]
