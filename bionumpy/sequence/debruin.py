from .kmers import get_kmers, KmerEncoding
from ..encodings import DNAEncoding
from ..encoded_array import as_encoded_array, EncodedArray
from typing import List, Dict


class DeBruijnGraph:
    def __init__(self, kmer_set, k):
        self._kmer_set = kmer_set
        print(self._kmer_set)
        self._kmer_encoding = KmerEncoding(DNAEncoding, k)
        self._k = k

    @classmethod
    def from_sequences(cls, sequences, k=31):
        kmers = get_kmers(as_encoded_array(sequences), k)
        return cls(set(kmers.ravel().raw()), k)

    def _get_previous(self, kmer):
        mask = (4**(self._k)-1)
        base = (kmer << 2) & mask
        return [base+i for i in range(4)]

    def _get_next(self, kmer):
        base = kmer >> 2
        return [base+(i << (2*(self._k-1))) for i in range(4)]

    def forward(self, kmer):
        kmer = as_encoded_array(kmer, self._kmer_encoding).raw()
        possible_next = self._get_next(kmer)
        return [self._kmer_encoding.to_string(n) for n in possible_next if n in self._kmer_set]

    def backward(self, kmer):
        kmer = as_encoded_array(kmer, self._kmer_encoding).raw()
        possible_next = self._get_previous(kmer)
        return [self._kmer_encoding.to_string(n) for n in possible_next if n in self._kmer_set]

    def get_contigs(self):
        pass

    def kmers(self):
        pass


class KallistoIndex:
    def __init__(self, lookup_table: Dict[int, List[int]]):
        self._lookup_table = lookup_table

    def from_debruin_graph(self, debruin_graph, sequences):
        pass

    def get_equivalence_classes(self, reads):
        pass
