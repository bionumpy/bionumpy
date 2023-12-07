import numpy as np
from typing import List, Union
from npstructures import RaggedArray, RaggedShape
from npstructures.util import unsafe_extend_right, unsafe_extend_left
from npstructures.raggedarray.raggedslice import ragged_slice

from bionumpy.encoded_array import EncodedArray, EncodedRaggedArray, change_encoding, as_encoded_array
from ..encodings.alphabet_encoding import DigitEncoding
from ..encodings import BaseEncoding
from ..string_array import StringArray


def int_to_str(number: int) -> str:
    number = np.asanyarray(number)
    L = np.log10(np.maximum(number, 1)).astype(int)+1
    digits = number // 10**np.arange(L)[::-1] % 10
    return EncodedArray(digits, DigitEncoding)


def _build_power_array(shape: RaggedShape, dots: np.ndarray = None) -> RaggedArray:
    """Build a ragged array where each row is the tenth power of the digit at that position

    If dots is specified, keep that space for a dot in the number

    Parameters
    ----------
    shape : RaggedShape
        The shape of the resulting string
    dots : np.ndarray
        Where the dots are located

    Returns
    -------
    RaggedArray
        The power of each digit in the corresponding string

    Examples
    --------
    FIXME: Add docs.

    """
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


def replace_inplace(number_text: EncodedArray, replace_from: str, replace_to: str) -> None:
    """ Replace inplace.

    Parameters
    ----------
    number_text : EncodedArray
    replace_from : str
    replace_to : str
    """
    number_text[number_text == replace_from] = replace_to


def str_to_int_with_missing(number_text: EncodedArray, missing_value=0) -> np.ndarray:
    return parse_with_missing(missing_value, number_text, str_to_int, int)

def str_to_float_with_missing(number_text: EncodedArray, missing_value=np.nan) -> np.ndarray:
    return parse_with_missing(missing_value, number_text, str_to_float, float)

def parse_with_missing(missing_value, number_text, parser, dtype):
    number_text = as_encoded_array(number_text)
    if np.all(number_text.lengths == 1):
        if np.all(number_text=='.'):
            return np.full(len(number_text), missing_value, dtype=dtype)
    mask = number_text.lengths > 0
    values = np.full(len(number_text), missing_value, dtype=dtype)
    values[mask] = parser(number_text[mask])
    return values


def str_to_int(number_text: EncodedArray, is_negative=None, is_positive=None) -> np.ndarray:
    """Convert strings in an EncodedRaggedArray to integers in numpy array

    Parameters
    ----------
    number_text : EncodedArray
        integer-strings in ragged array

    Returns
    -------
    np.ndarray
        integer values

    Examples
    --------
    FIXME: Add docs.

    """
    number_text = as_encoded_array(number_text).copy()
    if not isinstance(number_text, EncodedRaggedArray):

        number_text = as_encoded_array(number_text, target_encoding=DigitEncoding)
        powers = 10**np.arange(number_text.shape[-1])[::-1]
        return number_text.raw().dot(powers)
    else:
        number_text.ravel()
    if len(number_text) == 0:
        return np.array([], dtype=int)
    if is_negative is None:
        is_negative = number_text[:, 0] == "-"
        is_positive = number_text[:, 0] == "+"
    number_text[is_negative, 0] = "0"
    number_text[is_positive, 0] = "0"
    number_text = as_encoded_array(number_text, DigitEncoding)
    number_digits = RaggedArray(number_text.ravel().data, number_text._shape)
    powers = 10**_build_power_array(number_text._shape)
    signs = np.where(is_negative, -1, +1)
    return (number_digits*powers).sum(axis=-1)*signs


def _decimal_str_to_float(number_text: EncodedArray) -> np.ndarray:
    number_text = as_encoded_array(number_text)
    is_negative = number_text[:, 0] == "-"
    number_text[is_negative, 0] = "0" 
    dots = np.nonzero(number_text == ".")
    number_text[dots] = "0"
    number_text = as_encoded_array(number_text, DigitEncoding)
    power_array = _build_power_array(number_text._shape, dots=dots)
    powers = 10.**power_array
    number_digits = number_text.raw()
    base_numbers = (number_digits*powers).sum(axis=-1)
    row_indices, col_indices = dots
    exponents = np.zeros_like(number_text.lengths)
    exponents[row_indices] = number_text.lengths[row_indices] - col_indices-1
    powers = (10.**(exponents))
    signs = np.where(is_negative, -1, +1)
    return signs*base_numbers / powers


def _scientific_str_to_float(number_text: EncodedArray) -> np.ndarray:
    number_text = as_encoded_array(number_text)
    row, cols = np.nonzero(number_text == "e")
    decimal_text = ragged_slice(number_text, ends=cols)
    decimal_numbers = _decimal_str_to_float(decimal_text)
    power_text = ragged_slice(number_text, starts=cols+1)
    powers = str_to_int(power_text)
    return decimal_numbers*10.**powers


def str_to_float(number_text: EncodedRaggedArray) -> np.ndarray:
    """Convert strings representing floats to floats

    Parameters
    ----------
    number_text : EncodedRaggedArray
        The strings to be converted

    Returns
    -------
    np.ndarray
        Numpy array with the floats

    Examples
    --------
    FIXME: Add docs.

    """
    number_text = as_encoded_array(number_text)
    assert number_text.encoding == BaseEncoding
    # return np.array([float(row.to_string()) for row in number_text])

    scientific = np.any(number_text == "e", axis=-1)
    numbers = np.empty(len(number_text))
    if np.sum(scientific):
        numbers[scientific] = _scientific_str_to_float(number_text[scientific])
    if np.sum(~scientific):
        numbers[~scientific] = _decimal_str_to_float(number_text[~scientific])
    return numbers


def ints_to_strings(number: np.ndarray) -> EncodedRaggedArray:
    """Convert an array of ints into an ecoded ragged array holding their string representation

    Parameters
    ----------
    number : np.ndarray
        The numbers to be converted

    Returns
    -------
    EncodedRaggedArray
        The string representations in an EncodedRaggedArray

    Examples
    --------
    FIXME: Add docs.

    """
    number = np.asanyarray(number)
    is_negative = number < 0
    lengths = np.log10(np.maximum(np.abs(number), 1)).astype(int)+1
    shape = RaggedShape(lengths+is_negative)
    ragged_index = _build_power_array(shape)
    digits = np.abs(number)[:, np.newaxis] // 10**ragged_index % 10
    digits = EncodedRaggedArray(
        EncodedArray(digits.ravel(), DigitEncoding), digits._shape)
    digits = change_encoding(digits, BaseEncoding)
    #digits = as_encoded_array(digits, target_encoding=BaseEncoding)
    digits[is_negative, 0] = "-"
    return digits


def float_to_strings(floats: np.ndarray) -> EncodedRaggedArray:
    """Convert floats to strings

    This is actually very hard, need to write dragon4 or similar to
    numpy. For now use vanilla-speed

    Parameters
    ----------
    floats : np.ndarray
        floats to convert

    Returns
    -------
    EncodedRaggedArray
        strings in EncodedRaggedArray

    """
    return as_encoded_array([str(f) for f in floats])
    s = np.array2string(floats, max_line_width=10**15, threshold=10**15, separator=",").replace(" ", "")
    assert " " not in s, s
    b = np.frombuffer(bytes(s, encoding="ascii"), dtype=np.uint8)
    return split(EncodedArray(b[1:-1]), sep=",")


def int_lists_to_strings(int_lists: RaggedArray, sep: str = ",", keep_last: bool = False) -> EncodedRaggedArray:
    """Join the ints in each row into a string 

    Parameters
    ----------
    int_lists : RaggedArray
        RaggedArray of ints to be joined
    sep : str
        Separator character used to join the strings
    keep_last : bool
        Whether or not to keep the trailing separator on each line

    Returns
    -------
    EncodedRaggedArray
        EncodedRaggeadArray with the joined ints

    Examples
    --------
    FIXME: Add docs.

    """
    if len(sep) == 0:
        return EncodedRaggedArray(EncodedArray(int_lists.ravel(), DigitEncoding), int_lists._shape)
    int_strings = ints_to_strings(int_lists.ravel())
    lengths = RaggedArray(int_strings.lengths, int_lists._shape)
    joined = join(int_strings, sep=sep, keep_last=True)
    row_lens = lengths.sum(axis=-1)+int_lists.lengths
    ra = EncodedRaggedArray(joined, row_lens)
    if not keep_last:
        ra = ra[:, :-1]
    return ra


def join(sequences: EncodedRaggedArray, sep: str = "\t", keep_last: bool = False) -> EncodedArray:
    """Join a set of sequences in an EncodedRaggedArray into a single EncodedArray

    Parameters
    ----------
    sequences : EncodedRaggedArray
        A set of encoded sequences
    sep : str
        The character used to separate the sequences
    keep_last : bool
        Wheter or not to keep the trailing seperator

    Returns
    -------
    EncodedArray
        `sequences` joined by `sep`

    Examples
    --------

    """
    new_lengths = sequences.lengths+1
    new_array = sequences.__class__(
        EncodedArray(np.empty(shape=np.sum(new_lengths), dtype=np.uint8), sequences.encoding), new_lengths)
    new_array[:, :-1] = sequences
    new_array[:, -1] = sep
    if keep_last:
        return new_array.ravel()
    return new_array.ravel()[:-1]


def split(sequence: EncodedArray, sep: Union[str, List[str]] = ",") -> EncodedRaggedArray:
    """Split the sequence on `sep` character into a ragged array

    Parameters
    ----------
    sequence : EncodedArray
        The sequence to split
    sep : str
        The character to split on

    Returns
    -------
    EncodedRaggedArray
        EncodedRaggedArray where each row is a part after split

    Examples
    --------
    """
    us = unsafe_extend_right(sequence)
    if isinstance(sep, list):
        mask = (us == sep[0])
        for s in sep[1:]:
            mask |= (us == s)
    else:
        mask = us == sep
    mask[-1] = True
    sep_idx = np.flatnonzero(mask)
    lens = np.diff(unsafe_extend_left(sep_idx))
    lens[0] = sep_idx[0]+1
    ragged_array = EncodedRaggedArray(unsafe_extend_right(sequence), lens)
    return ragged_array[:, :-1]


def str_equal(sequences: EncodedRaggedArray, match_string: Union[str, EncodedRaggedArray]) -> np.ndarray:
    """Test if any of the sequences in `sequences` equals `match_string`

    Parameters
    ----------
    sequences : EncodedRaggedArray
        The set of sequences to test
    match_string : Union[str, EncodedRaggedArray]
        Either a single string to match or a set of sequences to match

    Returns
    -------
    np.ndarray
        Boolean array of which sequences matches the `match_string`

    """
    if isinstance(sequences, StringArray):
        return sequences == match_string

    sequences = as_encoded_array(sequences)

    if isinstance(sequences, EncodedRaggedArray) and isinstance(match_string, EncodedRaggedArray):
        return _str_equal_two_encoded_ragged_arrays(sequences, match_string)

    if isinstance(sequences, EncodedArray):
        return (len(sequences) == len(match_string)) and np.all(sequences==match_string)
    L = len(match_string)
    mask = (sequences.lengths == L)
    sequences.ravel()
    starts = sequences._shape.starts[mask]
    matrix = sequences.ravel()[starts[:, np.newaxis]+np.arange(L)]
    mask[mask] &= np.all(matrix == match_string, axis=-1)
    return mask


def _str_equal_two_encoded_ragged_arrays(sequences: EncodedRaggedArray, sequences_b: EncodedRaggedArray) -> np.ndarray:
    sequences = as_encoded_array(sequences)
    L = sequences_b.lengths
    mask = (sequences.lengths == L)
    mask[mask] &= (sequences[mask] == sequences_b[mask]).all(axis=-1)
    return mask
