from bionumpy.sequences import Sequence
from npstructures import RaggedArray
from bionumpy.encodings.alphabet_encoding import DigitArray
import numpy as np


def int_to_str(number):
    number = np.asanyarray(number)
    L = np.log10(number).astype(int)+1
    digits = number // 10**np.arange(L)[::-1] % 10
    return digits.view(DigitArray)


def ints_to_strings(number):
    number = np.asanyarray(number)
    lengths = np.log10(number).astype(int)+1
    total_lengths = np.sum(lengths)
    index_array = np.full(total_lengths, -1, dtype=int)
    index_array[np.cumsum(lengths)[:-1]] = lengths[1:]-1
    index_array[0] = lengths[0]-1
    np.cumsum(index_array, out=index_array)
    ragged_index = RaggedArray(index_array, lengths)
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
