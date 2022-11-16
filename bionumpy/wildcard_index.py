import numpy as np
from numpy.typing import ArrayLike
from collections import defaultdict
from bionumpy import EncodedRaggedArray


class WildCardIndex:
    def __init__(self, shape, letter_map):
        self._shape = shape
        self._letter_map = letter_map

    @classmethod
    def create_index(cls, sequences: EncodedRaggedArray) -> "WildCardIndex":
        shape = sequences.shape
        flat_sequences = sequences.ravel()
        letter_map = defaultdict(list)
        for index, letter in enumerate(flat_sequences):
            letter_map[str(letter)].append(index)
        letter_map = {key: np.array(value) for key, value in letter_map.items()}
        return cls(shape, letter_map)

    def get_pattern_indices(self, pattern: str) -> set:
        index_sets = []
        for index, letter in enumerate(pattern):
            if letter == ".":
                continue
            positions = self._letter_map[letter]
            positions -= index
            index_sets.append(set(positions))
        common_indices = set.intersection(*index_sets)
        sequence_indices = np.searchsorted(self._shape.starts, list(common_indices), side="right")-1
        print(np.unique(sequence_indices))
        return np.unique(sequence_indices)

    def _get_sequence_indices(self, flat_indices):
        sequence_indices = np.searchsorted(self._shape.starts, flat_indices)
        print(sequence_indices)