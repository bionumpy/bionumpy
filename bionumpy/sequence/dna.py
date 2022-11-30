from ..datatypes import Interval
from ..encoded_array import EncodedArray, EncodedRaggedArray
from ..encoded_array import as_encoded_array
from ..lookup import Lookup
import numpy as np


def _get_complement_lookup(alphabet_encoding):
    alphabet = alphabet_encoding.get_alphabet()
    complements = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}
    new_alphabet = "".join(complements[c] for c in alphabet)
    return Lookup(as_encoded_array(new_alphabet, alphabet_encoding),
                  alphabet_encoding)


def complement(_array):
    array = _array
    if isinstance(_array, EncodedRaggedArray):
        array = _array.ravel()
    assert isinstance(array, EncodedArray)
    lookup = _get_complement_lookup(array.encoding)
    new_data = lookup[array]
    if isinstance(_array, EncodedRaggedArray):
        new_data = EncodedRaggedArray(new_data, _array.shape)
    return new_data


def get_reverse_complement(array):
    return complement(array)[..., ::-1]


def get_strand_specific_sequences(encoded_array: EncodedArray, stranded_intervals: Interval):
    relevant_sequences = encoded_array[stranded_intervals.start:stranded_intervals.stop]
    rev_complimnet_seqs = get_reverse_complement(relevant_sequences)
    return np.where((stranded_intervals.strand.ravel() == "-")[:, np.newaxis],
                    rev_complimnet_seqs, relevant_sequences)
