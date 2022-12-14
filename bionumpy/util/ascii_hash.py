import numpy as np
from ..encoded_array import EncodedRaggedArray
from npstructures.mixin import NPSArray
from npstructures import RaggedArray, HashTable
n_letters = 128


def column_index_array(shape):
    row_lengths = shape[-1]
    size = np.sum(row_lengths)
    index_builder = np.ones(size + 1, dtype=int)
    index_builder[np.cumsum(row_lengths)] = 1-row_lengths
    index_builder[0] = 0
    np.cumsum(index_builder, out=index_builder)
    return RaggedArray(index_builder[:-1], shape)


def _get_power_array(n, mod):
    l = [1]
    for _ in range(n-1):
        l.append((l[-1]*n_letters % mod))
    return np.array(l).view(NPSArray)


def get_ascii_hash(encoded_array, mod):
    powers = _get_power_array(np.max(encoded_array.shape[-1]), mod)
    if isinstance(encoded_array, EncodedRaggedArray):
        col_indices = column_index_array(encoded_array.shape)
        powers = RaggedArray(powers[col_indices.ravel()], encoded_array.shape)

    return np.sum((powers*encoded_array.raw()) % mod, axis=-1) % mod


class AsciiHashTable:
    big_mod = 2**31-1
    def __init__(self, hash_table):
        self._hash_table = hash_table

    @classmethod
    def from_sequences(cls, encoded_ragged_array, modulo=103):
        ascii_hashes = get_ascii_hash(encoded_ragged_array, cls.big_mod)
        hash_table = HashTable(ascii_hashes, np.arange(len(encoded_ragged_array)), mod=modulo)
        return cls(hash_table)

    def __getitem__(self, encoded_array):
        hashes = get_ascii_hash(encoded_array, self.big_mod)
        return self._hash_table[hashes]
