import numpy as np
from numpy.typing import ArrayLike
from numbers import Number
from ..io.matrix_dump import Matrix
from ..util.typing import EncodedArrayLike
from ..encoded_array import EncodedArray


class EncodedCounts:
    alphabet: list
    counts: np.ndarray
    row_names: list = None

    def __init__(self, alphabet, counts, row_names=None):
        self.counts = counts
        self.alphabet = alphabet
        self.row_names = row_names

    def __str__(self):
        return "\n".join(f"{c}: {n}" for c, n in zip(self.alphabet, self.counts.T))

    def __repr__(self):
        return f'''EncodedCounts(alphabet={repr(self.alphabet)}, counts={repr(self.counts)}, row_names={repr(self.row_names)})'''

    def __eq__(self, other):
        if self.alphabet != other.alphabet:
            return False
        if not np.all(self.counts == other.counts):
            return False
        return True

    def __getitem__(self, idx):
        return self.counts[..., self.alphabet.index(idx)]

    def __add__(self, other):
        if isinstance(other, Number):
            o_counts = other
        else:
            assert self.alphabet==other.alphabet
            o_counts = other.counts
        return self.__class__(self.alphabet, self.counts+o_counts)

    def __radd__(self, other):
        if isinstance(other, Number):
            o_counts = other
        else:
            assert self.alphabet==other.alphabet
            o_counts = other.counts
        return self.__class__(self.alphabet, self.counts+o_counts)
    # return dataclasses.replace(self, counts=self.counts+o_counts)

    @property
    def proportions(self):
        s = self.counts.sum(axis=-1, keepdims=True)
        return np.where(s > 0, self.counts/s, 0)
    
    @property
    def proportion_matrix(self):
        s = self.counts.sum(axis=-1, keepdims=True)
        return Matrix(np.where(s > 0, self.counts/s, 0), col_names=self.alphabet)
            
        

    def get_count_for_label(self, label):
        return np.sum(self.counts[..., self.alphabet.index(l)] for l in label)

    @property
    def labels(self):
        return self.alphabet

    @classmethod
    def vstack(cls, counts):
        alphabet = counts[0].alphabet
        row_names = counts[0].row_names
        assert all(count.alphabet==alphabet for count in counts)
        ret = cls(alphabet, np.array([count.counts for count in counts], dtype="int"))
        if row_names is not None:
            ret.row_names = [count.row_names for count in counts]
        return ret

    def most_common(self, n=None):
        args = np.argsort(self.counts)[::-1]
        if n is not None:
            args = args[:n]
        return self.__class__(
            [self.alphabet[i] for i in args],
            self.counts[args])

    def as_dict(self):
        return dict(zip(self.alphabet, self.counts.T))


def count_encoded(values: EncodedArrayLike, weights: ArrayLike = None, axis: int = -1) -> EncodedCounts:
    """Count the occurances of encoded entries. Works on any encoding with finite alphabet

    Parameters
    ----------
    values : EncodedArrayLike
    weights : ArrayLike
    axis : int

    Returns
    -------
    EncodedCounts

    """
    weights2d = weights is not None and np.asanyarray(weights).ndim ==2
    if axis is None:
        values = values.ravel()
    if hasattr(values.encoding, "get_alphabet"):
        alphabet = values.encoding.get_alphabet()
    else:
        alphabet = values.encoding.get_labels()
    if isinstance(values, EncodedArray) and len(values.shape) == 1 and not weights2d:
        max_size = 1000000
        if len(values) > max_size and weights is None:
            counts = sum(np.bincount(values[i*max_size:(i+1)*max_size], minlength=len(alphabet))
                         for i in range(len(values)//max_size+1))
        else:
            counts = np.bincount(values, weights=weights, minlength=len(alphabet))
    elif axis == -1:
        if not weights2d:
            counts = np.array([np.bincount(row, weights=weights, minlength=len(alphabet)) for row in values])
        else:
            counts = np.array([np.bincount(values, weights=row, minlength=len(alphabet)) for row in weights])
    return EncodedCounts(alphabet, counts)
