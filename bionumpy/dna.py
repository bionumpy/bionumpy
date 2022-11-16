from .encoded_array import EncodedArray, as_encoded_array, EncodedRaggedArray
from .lookup import Lookup


def _get_complement_lookup(alphabet_encoding):
    alphabet = alphabet_encoding.get_alphabet()
    complements = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}
    new_alphabet = "".join(complements[c] for c in alphabet)
    return Lookup(as_encoded_array(new_alphabet, alphabet_encoding),
                  alphabet_encoding)


def complement(_array):
    array = _array
    if isinstance(_array, EncodedRaggedArray):
        array, _array.ravel()
    assert isinstance(array, EncodedArray)
    lookup = _get_complement_lookup(array.encoding)
    new_data = lookup[array]
    if isinstance(_array, EncodedRaggedArray):
        new_data = EncodedRaggedArray(new_data, _array.shape)
    return new_data


def reverse_compliment(array):
    return complement(array)[..., ::-1]
