import numpy as np
from npstructures import ragged_slice
from .strops import ints_to_strings, join, split, str_to_float
from ..encoded_array import as_encoded_array
from dataclasses import dataclass


@dataclass
class Matrix:
    row_names: list
    col_names: list
    data: np.ndarray


def parse_matrix(text, field_type=float, colname_type=str, rowname_type=str, sep='\t'):
    assert colname_type == str
    assert field_type == float
    text = as_encoded_array(text)
    line_endings = np.flatnonzero(text == '\n')
    if colname_type is not None:
        col_names = split(text[:line_endings[0]], sep)
        text = text[line_endings[0]+1:]
    else:
        col_names = None
    seps = np.flatnonzero((text == sep) | (text == '\n'))
    starts = np.insert(seps[:-1], 0, -1)+1
    ends = seps
    if rowname_type is None:
        numbers = str_to_float(ragged_slice(text, starts, ends))
        return Matrix(None, col_names, numbers.reshape(-1, len(col_names)))


def matrix_to_csv(matrix, header=None, sep=",", row_names=None):
    assert np.issubdtype(matrix.dtype,  np.integer)
    entries = as_encoded_array(ints_to_strings(matrix.ravel()))
    if header is not None:
        entries = np.concatenate((as_encoded_array(header), entries))

    lens = (entries.lengths+1).reshape(-1, matrix.shape[-1])
    line_endings = np.cumsum(lens.sum(axis=-1))
    joined = join(entries, sep, keep_last=True)
    joined[line_endings-1] = "\n"
    return joined
