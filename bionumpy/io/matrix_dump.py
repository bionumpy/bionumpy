import numpy as np
from npstructures import ragged_slice
from .strops import ints_to_strings, join, split, str_to_float, str_to_int
from ..encoded_array import as_encoded_array
from dataclasses import dataclass


@dataclass
class Matrix:
    data: np.ndarray
    row_names: list = None
    col_names: list = None


def read_matrix(filename, *args, **kwargs):
    return parse_matrix(open(filename).read(), *args, **kwargs)


def parse_matrix(text, field_type=float, colname_type=str, rowname_type=str, sep='\t'):
    assert colname_type == str
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
    n_cols = len(col_names)
    row_names = None
    if rowname_type is not None:
        row_names = ragged_slice(text, starts[::n_cols], ends[::n_cols])
        starts = starts.reshape(-1, n_cols)[:, 1:].ravel()
        ends = ends.reshape(-1, n_cols)[:, 1:].ravel()
        col_names = col_names[1:]
    f = str_to_int if field_type == int else str_to_float
    numbers = f(ragged_slice(text, starts, ends))
    return Matrix(numbers.reshape(-1, len(col_names)), row_names, col_names)


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
