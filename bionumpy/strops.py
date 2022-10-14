from bionumpy.sequences import Sequence, Sequences, as_encoded_sequence_array, ASCIIText
from npstructures import RaggedArray, RaggedShape
from bionumpy.encodings.alphabet_encoding import DigitArray
from npstructures.util import unsafe_extend_right, unsafe_extend_left
import numpy as np


def int_to_str(number):
    number = np.asanyarray(number)
    L = np.log10(number).astype(int)+1
    digits = number // 10**np.arange(L)[::-1] % 10
    return digits.view(DigitArray)



def _build_power_array(shape, dots=None):
    total_lengths = shape.ends[-1]
    lengths = shape.lengths
    index_array = np.full(total_lengths, -1, dtype=int)
    offset = 0 
    if dots is not None:
        index_array[dots] = 0
        offset = 1
    index_array[np.cumsum(lengths)[:-1]] = lengths[1:]-1-offset
    index_array[0] = lengths[0]-1-offset
    np.cumsum(index_array, out=index_array)
    return RaggedArray(index_array, shape)


def str_to_int(number_text):
    number_text = as_encoded_sequence_array(number_text, DigitArray)
    number_digits = RaggedArray(np.asarray(number_text.ravel()),
                                number_text.shape)
    powers = 10**_build_power_array(number_text.shape)
    return (number_digits*powers).sum(axis=-1)


def _decimal_str_to_float(number_text):
    number_text = as_encoded_sequence_array(number_text, ASCIIText)
    flat = number_text._data
    dots = np.flatnonzero(flat == ".")
    flat[dots] = "0"
    number_text = as_encoded_sequence_array(number_text, DigitArray)
    powers = 10**_build_power_array(number_text.shape, dots=dots)
    number_digits = RaggedArray(np.asarray(number_text.ravel()),
                                number_text.shape)
    base_numbers = (number_digits*powers).sum(axis=-1)
    _, col_indices =  number_text.shape.unravel_multi_index(dots)
    powers = (10**(number_text.shape.lengths-col_indices-1))
    return base_numbers / powers

def _scientific_str_to_float(number_text):
    pass

def str_to_float(number_text):
    return _decimal_str_to_float(number_text)
    pass

def ints_to_strings(number):
    number = np.asanyarray(number)
    lengths = np.log10(number).astype(int)+1
    shape = RaggedShape(lengths)
    ragged_index = _build_power_array(shape)
    digits = number[:, np.newaxis] // 10**ragged_index % 10
    return RaggedArray(digits.ravel().view(DigitArray), digits.shape)


def join(sequences, sep="\t", keep_last=False):
    new_lengths = sequences.shape.lengths+1
    new_array = sequences.__class__(np.empty(shape=np.sum(new_lengths), dtype=np.uint8).view(sequences.ravel().__class__), new_lengths)
    new_array[:, :-1] = sequences
    new_array[:, -1] = sep
    if keep_last:
        return new_array.ravel()
    return new_array.ravel()[:-1]


def split(sequence, sep=","):
    mask = unsafe_extend_right(sequence) == sep
    mask[-1] = True
    sep_idx = np.flatnonzero(mask)
    lens = np.diff(unsafe_extend_left(sep_idx))
    lens[0] = sep_idx[0]+1
    ragged_array = Sequences(unsafe_extend_right(sequence), lens)
    return ragged_array[:, :-1]


def str_equal(sequences, s):
    L = len(s)
    mask = (sequences.shape.lengths == L)
    starts = sequences.shape.starts[mask]
    matrix = sequences.ravel()[starts[:, np.newaxis]+np.arange(L)]
    mask[mask] &= np.all(matrix == s, axis=-1)
    return mask
