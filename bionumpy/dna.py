from bionumpy.sequences import EncodedArray, as_encoded_sequence_array, Sequences


def _get_complement_lookup(alphabet_encoding):
    alphabet = alphabet_encoding.get_alphabet()
    complements = {"a": "t", "g": "c", "c": "g", "t": "a"}
    new_alphabet = "".join(complements[c] for c in alphabet)
    return as_encoded_sequence_array(new_alphabet, alphabet_encoding)


def complement(_array):
    array = _array
    if isinstance(_array, Sequences):
        array, _array.ravel()
    assert isinstance(array, EncodedArray)
    lookup = _get_complement_lookup(array.encoding)
    new_data = lookup[array]
    if isinstance(_array, Sequences):
        new_data = Sequences(new_data, _array.shape)
    return new_data


def reverse_compliment(array):
    return complement(array)[..., ::-1]
