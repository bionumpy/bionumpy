import numpy as np
from numbers import Number
from .encodings.base_encoding import (BaseEncoding, Encoding,
                                      NumericEncoding, ASCIIEncoding)
from npstructures import RaggedArray


class EncodedArray(np.ndarray):
    """ 
    Class for data that could be written as characters, but is represented numpy arrays
    """
    encoding = None

    @classmethod
    def from_string(cls, s: str) -> "EncodedArray":
        return cls(shape=(len(s),), dtype=np.uint8, buffer=bytes(s, encoding="ascii"))

    def to_string(self) -> str:
        return "".join(chr(n) for n in self)

    def __str__(self) -> str:
        text = self.encoding.decode(self)
        if len(self.shape) == 0:
            return chr(int(self))
        if len(self.shape) == 1:
            return "".join(chr(n) for n in text[:20]) + "..."*(len(text)>20)
        a = np.array(["".join(chr(n) for n in seq[:20]) for seq in self.reshape(-1, self.shape[-1])]).reshape(self.shape[:-1])[:20]
        return str(a)

    def __view_scalar(self, elem):
        if isinstance(elem, Number):
            elem = np.array(elem)
        return elem.view(self.__class__)

    def __getitem__(self, idx):
        return self.__view_scalar(super().__getitem__(idx))

    def __setitem__(self, idx, value):
        value = as_encoded_sequence_array(value, self.__class__)
        super().__setitem__(idx, value)

    def __iter__(self):
        return (self.__view_scalar(elem) for elem in super().__iter__())

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method == "__call__" and ufunc in (np.equal, np.not_equal):
            return ufunc(*(np.asarray(as_encoded_sequence_array(a, self.__class__)) for a in inputs))
        return super().__array_ufunc__(ufunc, method, *inputs, **kwargs)

    def __array_function__(self, func, types, args, kwargs):
        if not all(issubclass(t, self.__class__) for t in types):
            return NotImplemented
        if func == np.bincount:
            return np.bincount(np.asarray(args[0]), *args[1:], **kwargs)
        if func == np.concatenate:
            return func([np.asarray(e) for e in args[0]]).view(self.__class__)
        elif func in (np.append, np.insert, np.lib.stride_tricks.sliding_window_view):
            return func(np.asarray(args[0]), *args[1:], **kwargs).view(self.__class__)
        
        return NotImplemented
        return super().__array_function__(func, types, args, kwargs)

    def ___array_finalize__(self, obj):
        if self.dtype == bool:
            return self
        if isinstance(self, obj.__class__):
            return self
        return self.view(self.__class__)

    def __array_wrap__(self, out_arr, context=None):
        if out_arr.dtype == bool:
            return np.asarray(out_arr, dtype=bool)
        return out_arr

    def __eq__(self, other):
        if isinstance(other, int):
            other = np.asanyarray(other)
        elif isinstance(other, (str, list)):
            other = as_encoded_sequence_array(other, self.__class__)
        return np.asarray(self) == np.asarray(other)


class ASCIIText(EncodedArray):
    encoding = BaseEncoding


class QualityEncoding:
    @classmethod
    def encode(cls, byte_array):
        return (byte_array-ord("!")).view(Quality)

    @classmethod
    def decode(cls, obj):
        return (obj+ord("!")).view(ASCIIText)


class Quality(EncodedArray):
    encoding = QualityEncoding

    def __str__(self):
        return "".join(chr(c) for c in self.decode())

    @classmethod
    def encode(cls, byte_array):
        return (byte_array-ord("!")).view(cls)

    def decode(self):
        return (self+ord("!")).view(ASCIIText)


class Sequences(RaggedArray):
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        assert isinstance(self._data, EncodedArray), self._data
        inputs = [as_encoded_sequence_array(i, self._data.__class__) for i in inputs]
        kwargs = {key: as_encoded_sequence_array(val, self._data.__class__) for key, val in kwargs.items()}
        ret = super().__array_ufunc__(ufunc, method, *inputs, **kwargs)
        if isinstance(ret._data, EncodedArray):
            return ret
        return RaggedArray(ret._data, ret.shape)


class NumericEncodedArray(EncodedArray):
    pass


class Sequence(ASCIIText):
    pass


def create_sequence_array_from_already_encoded_data(data: np.ndarray, encoding: Encoding) -> EncodedArray:
    if isinstance(data, RaggedArray):
        return Sequences(create_sequence_array_from_already_encoded_data(data.ravel(), encoding), data.shape)
    if isinstance(data, np.ndarray):
        return data.view(encoding)

    assert isinstance(data, (np.ndarray, RaggedArray))
    return Sequences(data._data, data.shape, encoding=encoding)


def as_numeric_encoded_array(data: np.ndarray, encoding: NumericEncoding):
    if isinstance(data, (str, list)):
        data = as_sequence_array(data)
    if isinstance(data, RaggedArray):
        return data.__class__(as_numeric_encoded_array(data.ravel(), encoding), data.shape)
    if isinstance(data, EncodedArray):
        assert data.encoding == BaseEncoding
        return encoding.encode(data)
    return np.asanyarray(data)


def isinstanceorsubclass(obj, types):
    return (isinstance(obj, type) and issubclass(obj, types)) or isinstance(obj, types)


def as_encoded_sequence_array(s, encoding: Encoding) -> EncodedArray:
    if isinstanceorsubclass(encoding, NumericEncoding):
        return as_numeric_encoded_array(s, encoding)
    s = as_sequence_array(s)
    if isinstance(s, RaggedArray):
        return s.__class__(as_encoded_sequence_array(s.ravel(), encoding), s.shape)
    if isinstance(encoding, type) and issubclass(encoding, EncodedArray):
        if isinstance(s, encoding):
            return s
        elif isinstance(s, ASCIIText):
            ret = encoding.encoding.encode(s)
            return ret.view(encoding)
        elif isinstance(encoding, type) and issubclass(encoding, ASCIIText):
            return s.encoding.decode(s).view(ASCIIText)
        else:
            assert False, (encoding, s.__class__, s.encoding)

    if s.encoding != encoding:
        e = encoding.encode(s)
        s = e.view(Sequence)
        s.encoding = encoding
        return s

    return s


def as_sequence_array(s) -> EncodedArray:#:, encoding=BaseEncoding):
    """
    Return a sequence array representation of s
    """
    if isinstance(s, EncodedArray):
        return s
    elif isinstance(s, np.ndarray):
        return s.view(ASCIIText)
    elif isinstance(s, RaggedArray):
        return Sequences(as_sequence_array(s._data), s.shape)
    elif isinstance(s, str):
        return Sequence.from_string(s)
    elif isinstance(s, list):
        return Sequences(as_sequence_array("".join(s)), [len(ss) for ss in s])
    else:
        raise Exception(f"Cannot convert {s} of class {type(s)} to sequence array")


def to_ascii(sequence_array, encoding=None):
    if isinstance(sequence_array, EncodedArray):
        assert encoding is None
        return sequence_array.encoding.decode(sequence_array).view(ASCIIText)
    if isinstance(sequence_array, np.ndarray):
        assert encoding is not None
        return encoding.decode(sequence_array).view(ASCIIText)
    elif isinstance(sequence_array, RaggedArray):
        return sequence_array.__class__(to_ascii(sequence_array.ravel(), encoding), sequence_array.shape)


def from_sequence_array(s) -> str:
    if isinstance(s, Sequence):
        return s.to_string()
    elif isinstance(s, RaggedArray):
        return [k.to_string() for k in s]

