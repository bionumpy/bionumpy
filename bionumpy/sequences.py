import numpy as np
from .encodings import BaseEncoding
from npstructures import RaggedArray


class Sequence(np.ndarray):
    encoding = BaseEncoding

    @classmethod
    def from_string(cls, s):
        return cls(shape=(len(s),), dtype=np.uint8, buffer=bytes(s, encoding="ascii"))


class Sequences(RaggedArray):
    def __init__(self, data, shape=None, encoding=BaseEncoding):
        super().__init__(data, shape, dtype=np.uint8)
        self.encoding = encoding

    @classmethod
    def from_strings(cls, sequences):
        return cls([Sequence.from_string(seq) for seq in sequences])

    def to_sequences(self):
        return ["".join(chr(i) for i in array) for array in self.tolist()]

    def __str__(self):
        strings = (
            "".join(chr(i) for i in array[:20]) + "..." * (len(array) > 20)
            for array in self.tolist()
        )
        seqs = ", ".join(seq for seq, _ in zip(strings, range(20)))
        trail = " ..." if len(self) > 20 else ""
        return f"Sequences({seqs}{trail})"

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        r = super().__array_ufunc__(ufunc, method, *inputs, **kwargs)
        if ufunc not in (np.equal, np.not_equal):
            return r

        return RaggedArray(r._data, r.shape)


def as_sequence_array(s, encoding=BaseEncoding):
    if isinstance(s, (Sequence, Sequences, np.ndarray, RaggedArray)):
        # assert s.encoding == encoding
        return s
    elif isinstance(s, str):
        return Sequence.from_string(s)
    elif isinstance(s, list):
        return Sequences.from_strings(s)
    else:
        raise Exception(f"Cannot convert {s} of class {type(s)} to sequence array")
