import numpy as np
from ..encoded_array import Encoding, as_encoded_array, EncodedArray, EncodedRaggedArray
from .exceptions import EncodingError
from ..util.ascii_hash import AsciiHashTable


class StringEncoding(Encoding):
    def __init__(self, sequences, modulo=None):
        self._seqeunces = as_encoded_array(sequences)
        self._modulo = modulo
        self._hash_table = AsciiHashTable.from_sequences(self._seqeunces)

    def get_labels(self):
        return self._seqeunces.tolist()

    def to_string(self, n):
        return self._seqeunces[int(n)].to_string()

    def encode(self, encoded_ragged_array):
        encoded_ragged_array = as_encoded_array(encoded_ragged_array)
        is_flat = isinstance(encoded_ragged_array, EncodedArray)
        if is_flat:
            encoded_ragged_array = EncodedRaggedArray(encoded_ragged_array, [len(encoded_ragged_array)])
        try:
            hashes = self._hash_table[encoded_ragged_array]
        except IndexError as e:
            raise EncodingError('String encoding failed') from e
        if is_flat:
            hashes = np.squeeze(hashes)
        return EncodedArray(hashes, self)

    def decode(self, encoded_array):
        if isinstance(encoded_array, EncodedArray):
            data = encoded_array.raw()
        else:
            data = encoded_array
        return self._seqeunces[data]

    def __repr__(self):
        return f'StringEncoding({self._seqeunces.tolist()})'

    def __eq__(self, other):
        if not isinstance(other, StringEncoding):
            return False
        my_shape, other_shape = (o._seqeunces.shape[-1] for o in (self, other))
        if len(my_shape) != len(other_shape):
            return False
        if np.any(my_shape != other_shape):
            return False
        return np.all(self._seqeunces == other._seqeunces) and self._modulo == other._modulo
