import numpy as np
import cupy as cp
from ..encodings import BaseEncoding
from npstructures import RaggedArray

#class CPSequence(cp.ndarray):
class CPSequence():
    encoding = BaseEncoding

    def __init__(self, shape, dtype, data):
        self._ndarray = cp.ndarray(shape=shape, dtype=dtype, memptr=data)

    def __getitem__(self, index):
        return self._ndarray[index]

    @classmethod
    def from_string(cls, s):
        return cls(shape=(len(s),), dtype=np.uint8, memptr=bytes(s, encoding="ascii"))

    def __repr__(self):
        self._ndarray = self.encoding.decode(self._ndarray)
        if len(self.shape) == 1:
            return "Sequence(" + "".join(chr(n) for n in self[:20]) + ")"
        return "Sequence(" + str(cp.array(["".join(chr(n) for n in seq[:20]) for seq in self._ndarray.reshape(-1, self._ndarray.shape[-1])]).reshape(self._ndarray.shape[:-1])[:20]) + ")"

    def __str__(self):
        if len(self._ndarray.shape) == 1:
            return "".join(chr(n) for n in self.encoding.decode(self._ndarray)[:10])
        return self.__repr__()

    @classmethod
    def from_array(cls, a):
        if isinstance(a, CPSequence):
            return cls(a._ndarray.shape, a._ndarray.dtype, a._ndarray.data)
        return cls(a.shape, a.dtype, a.data)


class CPSequences(RaggedArray):
    def __init__(self, data, shape=None, encoding=BaseEncoding):
        super().__init__(data, shape, dtype=np.uint8)
        self._data = CPSequence.from_array(self._data)
        self._data.encoding = encoding
        self.encoding = encoding

    def ravel(self):
        d = super().ravel()
        s = CPSequence.from_array(d)
        s.encoding = self.encoding
        return s

    @classmethod
    def from_strings(cls, sequences):
        return cls([CPSequence.from_string(seq) for seq in sequences])

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
