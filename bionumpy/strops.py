from .encoded_array import EncodedArray, as_encoded_array, EncodedRaggedArray
from npstructures import RaggedArray, RaggedShape
from bionumpy.encodings.alphabet_encoding import DigitEncoding
from npstructures.util import unsafe_extend_right, unsafe_extend_left
from npstructures.raggedarray.raggedslice import ragged_slice
import numpy as np


def int_to_str(number):
    number = np.asanyarray(number)
    L = np.log10(number).astype(int)+1
    digits = number // 10**np.arange(L)[::-1] % 10
    return EncodedArray(digits, DigitEncoding)

def _build_power_array(shape, dots=None):
    total_lengths = shape.ends[-1]
    lengths = shape.lengths
    index_array = np.full(total_lengths, -1, dtype=int)
    offset_0, offset_rest = (0, 0)
    if dots is not None:
        index_array[shape.ravel_multi_index(dots)] = 0
        offset = np.zeros(len(lengths), dtype=int)
        offset[dots[0]] = 1
        offset_0, offset_rest = (offset[0], offset[1:])
    index_array[np.cumsum(lengths)[:-1]] += lengths[1:]-offset_rest
    index_array[0] += lengths[0]-offset_0
    np.cumsum(index_array, out=index_array)
    return RaggedArray(index_array, shape)


def str_to_int(number_text):
    number_text = as_encoded_array(number_text)
    is_negative = number_text[:, 0] == "-"
    number_text[is_negative, 0] = "0"
    number_text = as_encoded_array(number_text, DigitEncoding)
    number_digits = RaggedArray(number_text.ravel().data, number_text.shape)
    powers = 10**_build_power_array(number_text.shape)
    signs = np.where(is_negative, -1, +1)
    return (number_digits*powers).sum(axis=-1)*signs


def _decimal_str_to_float(number_text):
    number_text = as_encoded_array(number_text)
    is_negative = number_text[:, 0] == "-"
    number_text[is_negative, 0] = "0" 
    dots = np.nonzero(number_text == ".")
    number_text[dots] = "0"
    number_text = as_encoded_array(number_text, DigitEncoding)
    power_array = _build_power_array(number_text.shape, dots=dots)
    powers = 10.**power_array
    number_digits = number_text.raw()
    base_numbers = (number_digits*powers).sum(axis=-1)
    row_indices, col_indices = dots# number_text.shape.unravel_multi_index(dots)
    exponents = np.zeros_like(number_text.shape.lengths)
    exponents[row_indices] = number_text.shape.lengths[row_indices] - col_indices-1
    powers = (10.**(exponents))
    signs = np.where(is_negative, -1, +1)
    return signs*base_numbers / powers

def _scientific_str_to_float(number_text):
    number_text = as_encoded_array(number_text)
    row, cols = np.nonzero(number_text == "e")
    decimal_text = ragged_slice(number_text, ends=cols)
    decimal_numbers = _decimal_str_to_float(decimal_text)
    power_text = ragged_slice(number_text, starts=cols+1)
    powers = str_to_int(power_text)
    return decimal_numbers*10.**powers


def str_to_float(number_text):
    number_text = as_encoded_array(number_text)
    scientific = np.any(number_text == "e", axis=-1)
    numbers = np.empty(len(number_text))
    if np.sum(scientific):
        numbers[scientific] = _scientific_str_to_float(number_text[scientific])
    if np.sum(~scientific):
        numbers[~scientific] = _decimal_str_to_float(number_text[~scientific])
    return numbers
# return _decimal_str_to_float(number_text)


def ints_to_strings(number):
    number = np.asanyarray(number)
    lengths = np.log10(number).astype(int)+1
    shape = RaggedShape(lengths)
    ragged_index = _build_power_array(shape)
    digits = number[:, np.newaxis] // 10**ragged_index % 10
    return EncodedRaggedArray(
        EncodedArray(digits.ravel(), DigitEncoding), digits.shape)


def join(sequences, sep="\t", keep_last=False):
    new_lengths = sequences.shape.lengths+1
    new_array = sequences.__class__(
        EncodedArray(np.empty(shape=np.sum(new_lengths), dtype=np.uint8), sequences.encoding), new_lengths)
    new_array[:, :-1] = sequences
    new_array[:, -1] = sep
    if keep_last:
        return new_array.ravel()
    return new_array.ravel()[:-1]


def split(sequence, sep=","):
    us = unsafe_extend_right(sequence)
    mask = us == sep
    mask[-1] = True
    sep_idx = np.flatnonzero(mask)
    lens = np.diff(unsafe_extend_left(sep_idx))
    lens[0] = sep_idx[0]+1
    ragged_array = EncodedRaggedArray(unsafe_extend_right(sequence), lens)
    return ragged_array[:, :-1]


def str_equal(sequences, s):
    L = len(s)
    mask = (sequences.shape.lengths == L)
    starts = sequences.shape.starts[mask]
    matrix = sequences.ravel()[starts[:, np.newaxis]+np.arange(L)]
    mask[mask] &= np.all(matrix == s, axis=-1)
    return mask
