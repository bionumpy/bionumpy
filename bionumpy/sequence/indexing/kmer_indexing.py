from ...encoded_array import EncodedRaggedArray, as_encoded_array
from .. import get_kmers
import numpy as np
from collections import defaultdict


class KmerIndex:
    def __init__(self, k, lookup, sequences_encoding):
        """

        """
        self._k = k
        self._lookup = lookup
        self._sequences_encoding = sequences_encoding

    def __repr__(self):
        return f"{self._k}-merIndex of sequences with {self._sequences_encoding}"

    @property
    def k(self):
        return self._k

    @classmethod
    def create_index(cls, sequences: EncodedRaggedArray, k: int) -> "KmerIndex":
        """
        Create a k-mer index of sequences with a given k

        Parameters
        ----------
        sequences: EncodedRaggedArray
        k: int

        Returns
        -------
        KmerIndex

        """
        kmers = get_kmers(sequences, k).raw()
        unique_kmers = np.unique(kmers.ravel())
        lookup = defaultdict(list)
        for kmer in unique_kmers:
            lookup[int(kmer)] = cls._get_index_for_kmer(kmers, kmer)
        return cls(k, lookup, sequences.encoding)

    @classmethod
    def _get_index_for_kmer(cls, kmers, kmer):
        return np.flatnonzero(np.any(kmers == kmer, axis=-1))

    def get_indices(self, kmer):
        if isinstance(kmer, str):
            assert len(kmer) == self._k
            encoded_kmer = get_kmers(sequence=as_encoded_array(kmer,self._sequences_encoding), k=self._k).raw()
            return self._lookup[int(encoded_kmer)]
        else:
            return self._lookup[int(kmer)]


class KmerLookup:

    index_class = KmerIndex

    def __init__(self, kmer_index: index_class, sequences: EncodedRaggedArray):
        self._kmer_index = kmer_index
        self._sequences = sequences

    def __repr__(self):
        return f"Lookup on {self._kmer_index.k}-merIndex of {len(self._sequences)} sequences"

    @classmethod
    def create_lookup(cls, sequences: EncodedRaggedArray, *args, **kwargs) -> "KmerLookup":
        index = cls.index_class.create_index(sequences=sequences, *args, **kwargs)
        return cls(index, sequences)

    def get_sequences(self, kmer):
        kmer_indices = self._kmer_index.get_indices(kmer)
        return self._sequences[kmer_indices]
