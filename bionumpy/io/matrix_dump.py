import numpy as np
from .strops import ints_to_strings, join
from ..encoded_array import as_encoded_array


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
