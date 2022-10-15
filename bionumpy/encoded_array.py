from npstructures import RaggedArray
from numbers import Number
from .encodings.base_encoding import BaseEncoding
from .encodings import Encoding, NumericEncoding
from .util import is_subclass_or_instance
import numpy as np


class EncodedRaggedArray(RaggedArray):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert isinstance(self._data, EncodedArray)

    @property
    def encoding(self):
        return self._data.encoding

    def raw(self):
        return RaggedArray(self._data.raw(), self.shape)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        assert isinstance(self._data, EncodedArray), self._data
        inputs = [as_encoded_array(i, self._data.encoding).raw() for i in inputs]
        kwargs = {key: as_encoded_array(val, self._data.encoding).raw() for key, val in kwargs.items()}
        print(ufunc, inputs, kwargs)
        ret = super().__array_ufunc__(ufunc, method, *inputs, **kwargs)
        if isinstance(ret._data, EncodedArray):
            return EncodedRaggedArray(ret._data, ret.shape)
        return ret


class EncodedArray(np.lib.mixins.NDArrayOperatorsMixin):
    """ 
    Class for data that could be written as characters, but is represented numpy arrays
    """

    encoding = None

    def __init__(self, data, encoding=BaseEncoding):
        if isinstance(data, EncodedArray):
            assert data.encoding == encoding
            data = data.data
        self.encoding = encoding
        self.data = np.asarray(data)

    def __len__(self):
        return len(self.data)

    def raw(self):
        return self.data

    def to_string(self):
        return "".join([chr(c) for c in self.encoding.decode(self.data)])

    def reshape(self, *args, **kwargs):
        return self.__class__(self.data.reshape(*args, **kwargs), self.encoding)

    @property
    def size(self):
        return self.data.size

    @property
    def shape(self):
        return self.data.shape

    @property
    def strides(self):
        return self.data.strides

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
        self.data.__setitem__(idx, value.data)

    def __iter__(self):
        return (self.__class__(element, self.encoding) for element in self.data)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method == "__call__" and ufunc in (np.equal, np.not_equal):
            return ufunc(*(as_encoded_array(a, self.encoding).raw() for a in inputs))
        return NotImplemented

    def __array_function__(self, func, types, args, kwargs):
        if func == np.bincount:
            return np.bincount(args[0].data, *args[1:], **kwargs)
        if func == np.concatenate:
            return self.__class__(func([e.data for e in args[0]]), self.encoding)
        if func == np.where:
            return self.__class__(func(args[0], args[1].data, args[2].data), encoding = self.encoding)
        elif func in (np.append, np.insert, np.lib.stride_tricks.sliding_window_view, np.lib.stride_tricks.as_strided):
            print(func, args[0], args[1:], kwargs, flush=True)
            return self.__class__(func(args[0].data, *args[1:], **kwargs), self.encoding)
        
        return NotImplemented
        return super().__array_function__(func, types, args, kwargs)

    def ravel(self):
        return self.__class__(self.data.ravel(), self.encoding)

    def as_strided(self, *args, **kwargs):
        assert isinstance(self.data, np.ndarray) and not np.issubdtype(self.data.dtype, np.object_)
        return self.__class__(np.lib.stride_tricks.as_strided(self.data, *args, **kwargs), self.encoding)


def as_encoded_array(s, target_encoding: Encoding = BaseEncoding) -> EncodedArray:
    if isinstance(s, str):
        s = EncodedArray([ord(c) for c in s], BaseEncoding)

    elif isinstance(s, list):
        s = EncodedRaggedArray(
            EncodedArray([ord(c) for ss in s for c in ss]),
            [len(ss) for ss in s])
    if isinstance(s, RaggedArray):
        data = as_encoded_array(s.ravel(), target_encoding)
        if isinstance(data, EncodedArray):
            return EncodedRaggedArray(data, s.shape)
        return RaggedArray(data, s.shape)
    if isinstance(s, np.ndarray):
        assert is_subclass_or_instance(target_encoding, NumericEncoding), s
        return s
    assert isinstance(s, EncodedArray), (s, repr(s), type(s))
    if s.encoding == target_encoding:
        return s
    elif s.encoding == BaseEncoding:
        encoded = target_encoding.encode(s.data)
        if is_subclass_or_instance(target_encoding, NumericEncoding):
            return encoded
        return EncodedArray(encoded, target_encoding)
    elif target_encoding == BaseEncoding:
        return EncodedArray(s.encoding.decode(s.data), BaseEncoding)
    assert False, (str(s.encoding), str(target_encoding))


def from_encoded_array(encoded_array: EncodedArray) -> str:
    if isinstance(encoded_array, EncodedRaggedArray):
        return [from_encoded_array(row) for row in encoded_array]
    else:
        return "".join(chr(c) for c in encoded_array.encoding.decode(encoded_array.data))
