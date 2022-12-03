import numpy as np
from numpy.typing import ArrayLike
from ...encoded_array import EncodedRaggedArray
from functools import reduce

from .kmer_indexing import KmerLookup


class WildCardIndex:
    def __init__(self, shape, letter_map):
        self._shape = shape
        self._letter_map = letter_map

    @classmethod
    def create_index(cls, sequences: EncodedRaggedArray) -> "WildCardIndex":
        shape = sequences._shape
        flat_sequences = sequences.ravel()
        letter_map = {letter: np.flatnonzero(flat_sequences == letter) for letter in sequences.encoding.get_labels()}
        return cls(shape, letter_map)

    def get_indices(self, pattern: str) -> ArrayLike:
        index_sets = (self._letter_map[letter] - index for index, letter in enumerate(pattern) if letter != ".")
        common_indices = reduce(np.intersect1d, index_sets)
        sequence_indices = np.searchsorted(self._shape.starts, list(common_indices), side="right")-1
        mask = common_indices + len(pattern) <= self._shape.ends[sequence_indices]
        return np.unique(sequence_indices[mask])


class WildCardLookup(KmerLookup):
    index_class = WildCardIndex

    def __repr__(self):
        return f"Lookup on WildcardIndex of {len(self._sequences)} sequences"



