import numpy as np
from .encodings.base_encoding import BaseEncoding, Encoding
from npstructures import RaggedArray


class EncodedArray(np.ndarray):
    """ 
    Class for data that could be written as characters, but is represented numpy arrays
    """
    
    encoding = BaseEncoding


    def decode(self):
        return self.encoding.decode(self)

    @classmethod
    def from_string(cls, s: str) -> "EncodedArray":
        return cls(shape=(len(s),), dtype=np.uint8, buffer=bytes(s, encoding="ascii"))

    def to_string(self) -> str:
        return "".join(chr(n) for n in self)

    def __str__(self) -> str:
        self = self.encoding.decode(self)
        if len(self.shape) == 1:
            return '"' + "".join(chr(n) for n in self[:20]) + '"'
        a = np.array(["".join(chr(n) for n in seq[:20]) for seq in self.reshape(-1, self.shape[-1])]).reshape(self.shape[:-1])[:20]
        return str(a)

    def __repr__(self) -> str:
        return f"Sequence({np.asarray(self)}, {self.encoding.__class__.__name__})"

    def __getitem__(self, key) -> "EncodedArray":
        r = super().__getitem__(key)
        if isinstance(r, (Sequence)):
            r.encoding = self.encoding

        return r

    def __array_function__(self, func, types, args, kwargs):
        if func == np.concatenate:
            return np.concatenate([np.asarray(e) for e in args[0]]).view(self.__class__)
        return super().__array_function__(func, types, args, kwargs)


    def __array_wrap__(self, out_arr, context=None):
        if out_arr.dtype == bool:
            return np.asarray(out_arr, dtype=bool)
        return out_arr

    def __eq__(self, other):
        if isinstance(other, int):
            other = np.asanyarray(other)
        else:
            other = as_encoded_sequence_array(other, self.encoding)
        return np.asarray(self) == np.asarray(other)


class ASCIIText(EncodedArray):
    @classmethod
    def encode(cls, byte_array):
        return byte_array

    def decode(self):
        return self


class Quality(EncodedArray):
    encoding = None
    def __str__(self):
        return "".join(chr(c) for c in self.decode())
    
    def __repr__(self):
        return f"{self.__class__.__name__}({repr(np.asarray(self))})"

    @classmethod
    def encode(cls, byte_array):
        return (byte_array-ord("!")).view(cls)

    def decode(self):
        return (self+ord("!")).view(ASCIIText)

#     def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
#         inputs = [np.asarray(as_encoded_sequence_array(i, self.encoding)) if isinstance(i, (str, list)) else np.asarray(i)
#                   for i in inputs]
#         ret = getattr(ufunc, method)(*inputs, **kwargs)
#         ret = ret.view(self.__class__)
#         ret.encoding = self.encoding
#         return ret
# 
#         ret = super().__array_ufunc__(ufunc, method, *inputs, **kwargs)
#         print(ret)
#         return ret


class NumericEncodedArray(EncodedArray):
    pass

class Sequence(EncodedArray):
    pass

def create_sequence_array_from_already_encoded_data(data: np.ndarray, encoding: Encoding) -> EncodedArray:
    if isinstance(data, np.ndarray):
        return data.view(EncodedArray)

    assert isinstance(data, (np.ndarray, RaggedArray))
    return Seqeunces(data._data, data.shape, encoding=encoding)

def as_encoded_sequence_array(s, encoding: Encoding) -> EncodedArray:
    print("<<<<<<<<<<<<<<", s, type(s))
    s = as_sequence_array(s)

    if isinstance(s, RaggedArray):
        print("!!!!!!!!!!!!!!!!!!!!!", repr(s._data))
        return s.__class__(as_encoded_sequence_array(s.ravel(), encoding), s.shape)
    print("AS ENCODED SEQUNCE ARRAY", s.encoding, encoding)
    print(repr(s))
    if isinstance(encoding, type) and issubclass(encoding, EncodedArray):
        if isinstance(s, encoding):
            return s
        else:
            ret =  encoding.encode(s)
            print("+++", repr(ret))
            return ret
                  

    if s.encoding != encoding:
        print(s.encoding, encoding)
        e = encoding.encode(s)
        s = e.view(Sequence)
        s.encoding = encoding 
        return s

    return s

def as_sequence_array(s) -> EncodedArray:#:, encoding=BaseEncoding):
    """
    Return a sequence array representation of s
    """
    print("??", repr(s), type(s))
    if isinstance(s, EncodedArray):
        return s
    elif isinstance(s, np.ndarray):
        return s.view(EncodedArray)
    elif isinstance(s, RaggedArray):
        return RaggedArray(as_sequence_array(s._data), s.shape)
    elif isinstance(s, str):
        print(s, Sequence.from_string(s))
        return Sequence.from_string(s)
    elif isinstance(s, list):
        return RaggedArray(as_sequence_array("".join(s)), [len(ss) for ss in s])
    else:
        raise Exception(f"Cannot convert {s} of class {type(s)} to sequence array")

def to_ascii(sequence_array):
    if isinstance(sequence_array, EncodedArray):
        return sequence_array.decode()
    elif isinstance(sequence_array, RaggedArray):
        return sequence_array.__class__(to_ascii(sequence_array.ravel()), sequence_array.shape)
                                        

def from_sequence_array(s) -> str:
    if isinstance(s, Sequence):
        return s.to_string()
    elif isinstance(s, RaggedArray):
        return [k.to_string() for k in s]
