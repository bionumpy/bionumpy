import numpy as np

from npstructures import ragged_slice


def find_sub_delimiters(data, starts, ends, delimiter=','):
    indices = np.flatnonzero(data == delimiter)
    indice_starts = np.searchsorted(indices, starts)
    indice_ends = np.searchsorted(indices, ends)
    return ragged_slice(indices, indice_starts, indice_ends)