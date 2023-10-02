from functools import reduce

import numpy as np

'''
P(X | M) = 
'''
'''
8 bits per kmer
cache optimized bloom filter

'''

def hash_function(offset):
    def f(kmer):
        return kmer ^ offset

    return f


class BloomFilter:
    def __init__(self, mask_size, hash_functions):
        self._hash_functions = hash_functions
        self._mask = np.zeros(mask_size, dtype=bool)

    @classmethod
    def from_m_and_k(cls, m, k, seed=12345):
        offsets = np.random.RandomState(seed).randint(0, m, k)
        return cls(m, [hash_function(offset, m) for offset in offsets])

    @classmethod
    def from_hash_functions_and_seqeuences(cls, hash_functions, sequence, mask_size):
        bloom_filter = cls(mask_size, hash_functions)
        bloom_filter.insert(sequence)
        return bloom_filter

    def insert(self, sequences):
        for function in self._hash_functions:
            self._mask[function(sequences) % self._mask.size] = True

    def __getitem__(self, idx):
        return reduce(np.logical_and, (self._mask[h(idx) % self._mask.size] for h in self._hash_functions))


class InterleavedBloomFilter:
    def __init__(self, hash_functions, mask):
        self._hash_functions = hash_functions
        self._mask = mask

    @classmethod
    def from_hash_functions_and_seqeuences(cls, hash_functions, sequences, mask_size):
        mask = np.zeros((mask_size, len(sequences)), dtype=bool)
        for function in hash_functions:
            for i, seqeunce in enumerate(sequences):
                mask[function(sequences), i] = True
        return cls(hash_functions, mask)

    def __getitem__(self, idx):
        kmer, seq_idx = idx
        return np.all([self._mask[h(kmer), seq_idx] for h in self._hash_functions],
                      axis=0)
