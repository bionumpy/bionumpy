import numpy as np
from .encodings import BaseEncoding
from npstructures import RaggedArray


class EncodedArray(np.ndarray):
    encoding = BaseEncoding

    @classmethod
    def from_string(cls, s):
        return cls(shape=(len(s),), dtype=np.uint8, buffer=bytes(s, encoding="ascii"))

    def __str__(self):
        self = self.encoding.decode(self)
        if len(self.shape) == 1:
            return "Sequence(" + "".join(chr(n) for n in self[:20]) + ")"
        return "Sequence(" + str(np.array(["".join(chr(n) for n in seq[:20]) for seq in self.reshape(-1, self.shape[-1])]).reshape(self.shape[:-1])[:20]) + ")"

    def __repr__(self):
        return f"Sequence({np.asarray(self)}, {self.encoding.__class__.__name__})"

    def __getitem__(self, key):
        r = super().__getitem__(key)
        if isinstance(r, (Sequence, Sequences)):
            r.encoding = self.encoding

        return r

    def __array_wrap__(self, out_arr, context=None):
        if out_arr.dtype == bool:
            return np.asarray(out_arr, dtype=bool)
        return out_arr


class NumericEncodedArray(EncodedArray):
    pass


class Sequence(EncodedArray):
    pass

class EncodedRaggedArray(RaggedArray):

    def __init__(self, data, shape=None, encoding=BaseEncoding):
        super().__init__(data, shape)
        self._data = self._data.view(Sequence)
        self._data.encoding = encoding
        self.encoding = encoding

    def ravel(self):
        d = super().ravel()
        s = d.view(Sequence)
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

    def __repr__(self):
        return f"Sequences({self._data}, self.shape, self.encoding)"

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        r = super().__array_ufunc__(ufunc, method, *inputs, **kwargs)
        if ufunc not in (np.equal, np.not_equal, np.less, np.less_equal, np.greater, np.greater_equal):
            return r

        return RaggedArray(r._data, r.shape)

    def __array_function__(self, func, types, args, kwargs):
        r = super().__array_function__(func, types, args, kwargs)
        if isinstance(r, (Sequence, Sequences)):
            r.encoding = self.encoding

        return r

    def __getitem__(self, key):
        r = super().__getitem__(key)
        if isinstance(r, (Sequence, Sequences)):
            r.encoding = self.encoding

        return r

class Sequences(EncodedRaggedArray):
    pass

class NumericEncodedRaggedArray(EncodedRaggedArray):
    pass
    

def create_sequence_array_from_already_encoded_data(data, encoding):
    assert isinstance(data, (np.ndarray, RaggedArray))
    return Seqeunces(data._data, data.shape, encoding=encoding)

def as_encoded_sequence_array(s, encoding):
    s = as_sequence_array(s)
    if s.encoding != encoding:
        assert s.encoding == BaseEncoding
        if isinstance(s, Sequences):
            return Sequences(encoding.encode(s.ravel()), s.shape, encoding=encoding)
        else:
            e = encoding.encode(s)
            s = e.view(Sequence)
            s.encoding = encoding 
            return s

    return s

def as_sequence_array(s):#:, encoding=BaseEncoding):
    """
    Return a sequence array representation of s
    """
    if isinstance(s, (Sequences, Sequence)):
        return s
    elif isinstance(s, np.ndarray):
        return s.view(Sequence)
    elif isinstance(s, RaggedArray):
        return Sequences(s._data, s.shape)
    elif isinstance(s, str):
        return Sequence.from_string(s)
    elif isinstance(s, list):
        return Sequences.from_strings(s)
    else:
        raise Exception(f"Cannot convert {s} of class {type(s)} to sequence array")
