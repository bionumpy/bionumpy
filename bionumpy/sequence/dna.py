from ..datatypes import Interval
from ..streams import streamable
from ..encodings import BaseEncoding, AlphabetEncoding
from ..encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from ..util.typing import EncodedArrayLike
from ..bnpdataclass.bnpdataclassfunction import apply_to_npdataclass
from .lookup import Lookup
import numpy as np

_complements = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}


def _get_complement_lookup(encoding):
    if isinstance(encoding, AlphabetEncoding):
        return _get_alphabet_encoding_complement_lookup(encoding)
    
    elif encoding == BaseEncoding:
        return _get_ascii_complement_lookup()
    raise ValueError(f"Invalid encoding for dna-complement: {encoding}")


def _get_alphabet_encoding_complement_lookup(alphabet_encoding):
    alphabet = alphabet_encoding.get_alphabet()
    new_alphabet = "".join(_complements[c] for c in alphabet)
    return Lookup(as_encoded_array(new_alphabet, alphabet_encoding),
                  alphabet_encoding)


def _get_ascii_complement_lookup():
    values = np.zeros(128, dtype=np.uint8)
    for key, value in _complements.items():
        values[ord(key)] = ord(value)
    return Lookup(EncodedArray(values, BaseEncoding))




def complement(_array):
    array = _array
    if isinstance(_array, EncodedRaggedArray):
        array = _array.ravel()
    assert isinstance(array, EncodedArray)
    lookup = _get_complement_lookup(array.encoding)
    new_data = lookup[array]
    if isinstance(_array, EncodedRaggedArray):
        new_data = EncodedRaggedArray(new_data, _array._shape)
    return new_data


@streamable()
@apply_to_npdataclass("sequence")
def get_reverse_complement(sequence: EncodedArrayLike) -> EncodedArrayLike:
    """Get the reverse complement of one or more DNA-sequences

    Parameters
    ----------
    sequence : EncodedArrayLike
        Dna sequence

    Returns
    -------
    EncodedArrayLike
        Reverse Complement

    """
    sequence = as_encoded_array(sequence)
    return complement(sequence)[..., ::-1]


@streamable()
def get_strand_specific_sequences(encoded_array: EncodedArray, stranded_intervals: Interval):
    relevant_sequences = encoded_array[stranded_intervals.start:stranded_intervals.stop]
    rev_complimnet_seqs = get_reverse_complement(relevant_sequences)
    return np.where((stranded_intervals.strand.ravel() == "-")[:, np.newaxis],
                    rev_complimnet_seqs, relevant_sequences)


@streamable()
def get_sequences(sequence: EncodedArray, intervals: Interval):
    return sequence[intervals.start:intervals.stop]
