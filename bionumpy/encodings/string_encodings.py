from ..encoded_array import Encoding, as_encoded_array, EncodedArray, EncodedRaggedArray
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
        hashes = self._hash_table[encoded_ragged_array]
        return EncodedArray(hashes, self)

    def decode(self, encoded_array):
        if isinstance(encoded_array, EncodedArray):
            data = encoded_array.raw()
        else:
            data = encoded_array
        return self._seqeunces[data]
