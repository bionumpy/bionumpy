from npstructures import RaggedArray
from numbers import Number
from .encodings.base_encoding import BaseEncoding
from .encodings import Encoding
import numpy as np


class EncodedRaggedArray(RaggedArray):
    @property
    def encoding(self):
        return self._data.encoding


class EncodedArray(np.lib.mixins.NDArrayOperatorsMixin):
    """ 
    Class for data that could be written as characters, but is represented numpy arrays
    """

    encoding = None

    def __init__(self, data, encoding=BaseEncoding):
        self.encoding = encoding
        self.data = np.asarray(data)

    def __len__(self):
        return len(self.data)

    @property
    def size(self):
        return self.data.size

    @property
    def dtype(self):
        return self.data.dtype

    def __str__(self) -> str:
        text = self.encoding.decode(self.data)
        if len(self.data.shape) == 0:
            return chr(int(self.data))
        if len(self.data.shape) == 1:
            return "".join(chr(n) for n in text[:20]) + "..."*(len(text)>20)
        a = np.array([str(self.__class__(seq, self.encoding)) for seq in self.data.reshape(-1, self.data.shape[-1])]).reshape(self.data.shape[:-1])[:20]
        return str(a)

    def __getitem__(self, idx):
        return self.__class__(self.data.__getitem__(idx), self.encoding)

    def __setitem__(self, idx, value):
        value = as_encoded_array(value, self.encoding)
        super().__setitem__(idx, value)

    def __iter__(self):
        return (self.__class__(element, self.encoding) for element in self.data)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method == "__call__" and ufunc in (np.equal, np.not_equal):
            return ufunc(*(as_encoded_array(a, self.encoding).data for a in inputs))
        return NotImplemented

    def __array_function__(self, func, types, args, kwargs):
        if not all(issubclass(t, self.__class__) for t in types):
            return NotImplemented
        if func == np.bincount:
            return np.bincount(args[0].data, *args[1:], **kwargs)
        if func == np.concatenate:
            return self.__class__(func([e.data for e in args[0]]), self.encoding)
        elif func in (np.append, np.insert, np.lib.stride_tricks.sliding_window_view):
            return self.__class__(func(args[0].data, *args[1:], **kwargs), self.encoding)
        
        return NotImplemented
        return super().__array_function__(func, types, args, kwargs)

    def ravel(self):
        return self.__class__(self.data.ravel(), self.encoding)


def as_encoded_array(s, target_encoding: Encoding = BaseEncoding) -> EncodedArray:
    if isinstance(s, str):
        s = EncodedArray([ord(c) for c in s], BaseEncoding)

    elif isinstance(s, list):
        print(s)
        s = EncodedRaggedArray(
            EncodedArray([ord(c) for ss in s for c in ss]),
            [len(ss) for ss in s])
    if isinstance(s, RaggedArray):
        return s.__class__(as_encoded_array(s.ravel(), target_encoding), s.shape)
    if s.encoding == target_encoding:
        return s
    elif s.encoding == BaseEncoding:
        return EncodedArray(target_encoding.encode(s.data), target_encoding)
    assert False, (s.encoding, target_encoding)
