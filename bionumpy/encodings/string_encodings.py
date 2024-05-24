from typing import Optional, List, Union

import numpy as np
from ..encoded_array import Encoding, as_encoded_array, EncodedArray, EncodedRaggedArray, encoded_array_from_nparray
from .exceptions import EncodingError
from ..string_array import StringArray
from ..util.ascii_hash import AsciiHashTable


class StringEncoding(Encoding):
    """
    Encodes strings into a numeric value corresponding to the index of the string in the input list.
    """
    def __init__(self, sequences: List[str], modulo: Optional[int] = None):
        self._seqeunces = as_encoded_array(sequences)
        self._modulo = modulo
        self._hash_table = AsciiHashTable.from_sequences(self._seqeunces)

    def get_labels(self) -> List[str]:
        return self._seqeunces.tolist()

    def to_string(self, n: int) -> str:
        """
        Get the string corresponding to the index n.

        Parameters
        ----------
        n: int

        Returns
        -------
        str
            The string corresponding to the index n.

        """
        return self._seqeunces[int(n)].to_string()

    def encode(self, encoded_ragged_array: Union[EncodedRaggedArray, StringArray, List[str]]) -> Union[EncodedArray, EncodedRaggedArray]:
        """
        Encode a string, list of strings or EncodedRaggedArray into an EncodedArray or EncodedRaggedArray of hashed strings.

        Parameters
        ----------
        encoded_ragged_array: Union[EncodedRaggedArray, StringArray, list[str]]
            The input data to encode.

        Returns
        -------
        Union[EncodedArray, EncodedRaggedArray]
            The encoded data.

        """
        if isinstance(encoded_ragged_array, StringArray):
            pass # encoded_ragged_array = encoded_array_from_nparray(encoded_ragged_array)
        else:
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

    def decode(self, encoded_array: Union[EncodedArray, np.ndarray]) -> Union[str, List[str]]:
        """
        Decode an EncodedArray or np.ndarray into a string or list of strings.

        Parameters
        ----------
        encoded_array

        Returns
        -------

        """

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
