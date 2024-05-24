import numpy as np
from ..encoded_array import EncodedRaggedArray
from npstructures.mixin import NPSArray
from npstructures import RaggedArray, HashTable

from ..string_array import StringArray

n_letters = 129


def column_index_array(shape):
    row_lengths = shape[-1]
    size = np.sum(row_lengths)
    index_builder = np.ones(size + 1, dtype=int)
    index_builder[np.cumsum(row_lengths)] = 1-row_lengths
    index_builder[0] = 0
    np.cumsum(index_builder, out=index_builder)
    return RaggedArray(index_builder[:-1], shape)


def _get_power_array(n, mod):
    '''
    (1*mod + k) * 128  =mod 128
    '''
    l = [1]
    for _ in range(n-1):
        l.append(((l[-1]*n_letters) % mod))
    return np.array(l).view(NPSArray)


def get_ascii_hash(encoded_array, mod):
    if len(encoded_array) == 0:
        return np.array([], dtype=int)
    if isinstance(encoded_array, StringArray):
        return get_ascii_hash_string_array(encoded_array, mod)
    return get_ascii_hash_ragged(encoded_array, mod)


def get_ascii_hash_ragged(encoded_array, mod):
    powers = _get_power_array(np.max(encoded_array.shape[-1]), mod)
    if isinstance(encoded_array, EncodedRaggedArray):
        col_indices = column_index_array(encoded_array.shape)
        powers = RaggedArray(powers[col_indices.ravel()], encoded_array.shape)
    elif isinstance(encoded_array, StringArray):

        powers = RaggedArray(powers, encoded_array.shape)
    return np.sum((powers * encoded_array.raw()) % mod, axis=-1) % mod


def get_ascii_hash_string_array(encoded_array, mod):
    array = encoded_array.as_bytes()
    x = _get_power_array(array.shape[-1], mod)
    A = (array * x % mod)
    # A *= x
    return np.sum(A, axis=-1) % mod
    # return (A @ x).ravel() % mod


class AsciiHashTable:
    big_mod = (2**31)-1

    def __init__(self, hash_table, sequences):
        self._hash_table = hash_table
        self._seqeunces = sequences

    @classmethod
    def from_sequences(cls, encoded_ragged_array, modulo=103):
        from collections import Counter
        ascii_hashes = get_ascii_hash(encoded_ragged_array, cls.big_mod)
        assert len(set(ascii_hashes)) == len(ascii_hashes), (len(set(ascii_hashes)), len(ascii_hashes))
        hash_table = HashTable(ascii_hashes, np.arange(len(encoded_ragged_array)), mod=modulo)
        return cls(hash_table, encoded_ragged_array)

    def __getitem__(self, encoded_array):
        hashes = get_ascii_hash(encoded_array, self.big_mod)
        try:
            values = self._hash_table[hashes]
        except IndexError as e:
            missing = ~self._hash_table.contains(hashes)
            raise IndexError(f'Error Looking for:\n {encoded_array[missing]}\nAvailable keys:\n{self._seqeunces}')
        return values
