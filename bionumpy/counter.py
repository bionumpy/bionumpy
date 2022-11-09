import numpy as np
from numbers import Number
import dataclasses

@dataclasses.dataclass
class EncodedCounts:
    alphabet: list
    counts: np.ndarray
    row_names: list = None

    def __str__(self):
        return "\n".join(f"{c}: {n}" for c, n in zip(self.alphabet, self.counts))

    def __getitem__(self, idx):
        return self.counts[..., self.alphabet.index(idx)]

    def __add__(self, other):
        if isinstance(other, Number):
            o_counts = other
        else:
            assert self.alphabet==other.alphabet
            o_counts = other.counts
        return dataclasses.replace(self, counts=self.counts+o_counts)

    def __radd__(self, other):
        if isinstance(other, Number):
            o_counts = other
        else:
            assert self.alphabet==other.alphabet
            o_counts = other.counts
        return dataclasses.replace(self, counts=self.counts+o_counts)

    def get_count_for_label(self, label):
        return np.sum(self.counts[..., self.alphabet.index(l)] for l in label)

    @classmethod
    def vstack(cls, counts):
        alphabet = counts[0].alphabet
        row_names = counts[0].row_names
        assert all(count.alphabet==alphabet for count in counts)
        ret =  cls(alphabet, np.array([count.counts for count in counts], dtype="int"))
        if row_names is not None:
            ret.row_names = [count.row_names for count in counts]
        return ret
            


def count_encoded(values, weights=None, axis=-1):
    if hasattr(values.encoding, "get_alphabet"):
        alphabet = values.encoding.get_alphabet()
    else:
        alphabet = values.encoding.get_labels()
    counts = np.bincount(values, weights=weights, minlength=len(alphabet))
    return EncodedCounts(alphabet, counts)
