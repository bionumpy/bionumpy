import numpy as np
from .encodings import BaseEncoding
from npstructures import RaggedArray


class Sequence(np.ndarray):
    encoding = BaseEncoding

    @classmethod
    def from_string(cls, s):
        return cls(shape=(len(s),), dtype=np.uint8, buffer=bytes(s, encoding="ascii"))

    def __repr__(self):
        self = self.encoding.decode(self)
        if len(self.shape) == 1:
            return "Sequence(" + "".join(chr(n) for n in self[:20]) + ")"
        return "Sequence(" + str(np.array(["".join(chr(n) for n in seq[:20]) for seq in self.reshape(-1, self.shape[-1])]).reshape(self.shape[:-1])[:20]) + ")"

    def __str__(self):
        if len(self.shape) == 1:
            return "".join(chr(n) for n in self.encoding.decode(self)[:10])
        return self.__repr__()

    @classmethod
    def from_array(cls, a):
        return cls(a.shape, a.dtype, a.data)


class Sequences(RaggedArray):
    def __init__(self, data, shape=None, encoding=BaseEncoding):
        super().__init__(data, shape, dtype=np.uint8)
        self._data = Sequence.from_array(self._data)
        self._data.encoding = encoding
        self.encoding = encoding

    def ravel(self):
        d = super().ravel()
        s = Sequence.from_array(d)
        s.encoding = self.encoding
        return s

    @classmethod
    def from_strings(cls, sequences):
        return cls([Sequence.from_string(seq) for seq in sequences])

    def to_sequences(self):
        return ["".join(chr(i) for i in array) for array in self.tolist()]

    def __str__(self):
        self = self.__class__(self.encoding.decode(self.ravel()),
                              self.shape)
        strings = (
            "".join(chr(i) for i in array[:20]) + "..." * (len(array) > 20)
            for array in self.tolist()
        )
        seqs = ", ".join(seq for seq, _ in zip(strings, range(20)))
        trail = " ..." if len(self) > 20 else ""
        return f"Sequences({seqs}{trail})"

    __repr__ = __str__

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        r = super().__array_ufunc__(ufunc, method, *inputs, **kwargs)
        if ufunc not in (np.equal, np.not_equal):
            return r

        return RaggedArray(r._data, r.shape)


def as_sequence_array(s, encoding=BaseEncoding):
    if isinstance(s, (Sequences, Sequence)):
        if s.encoding != encoding:
            assert s.encoding == BaseEncoding
            if isinstance(s, Sequences):
                return Sequences(encoding.encode(s.ravel()), s.shape)
            else:
                e = encoding.encode(s)
                s = Sequence.from_array(e)
                s.encoding=encoding
                return s
        return s
    elif isinstance(s, np.ndarray):
        s = encoding.encode(s)
        s = Sequence.from_array(s)
        s.encoding = encoding
        return s
    elif isinstance(s, RaggedArray):
        return s
    elif isinstance(s, str):
        return encoding.encode(Sequence.from_string(s))
    elif isinstance(s, list):
        s = Sequences.from_strings(s)
        return Sequences(encoding.encode(s.ravel()), s.shape, encoding=encoding)
    else:
        raise Exception(f"Cannot convert {s} of class {type(s)} to sequence array")
