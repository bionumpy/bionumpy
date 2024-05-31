"""
NB: This file contains EncodedArray, EncodedRaggedArray, as_encoded_array and BaseEncoding/Encoding classes.
For now, all these are depending on each other and needs to be in the same file in order to
avoid circular imports.
Should be refactored later.
"""
from numbers import Number

from npstructures import RaggedArray
from npstructures.mixin import NPSArray
from typing import Tuple, List, Union
import numpy as np
from abc import abstractmethod


class Encoding:
    @abstractmethod
    def encode(self, *args, **kwargs):
        return NotImplemented

    @abstractmethod
    def get_labels(self):
        pass

    def __call__(self, *args, **kwargs):
        return self.encode(*args, **kwargs)

    def is_base_encoding(self):
        return False

    def is_one_to_one_encoding(self):
        return False

    def is_numeric(self):
        return False


class OneToOneEncoding(Encoding):
    """Represents encodings that are one-to-one, i.e. where each element is
    encoded to one element and vice versa. This class is meant to be subclassed
    when implementing specific encodings."""

    def encode(self, data):
        assert hasattr(self, "_encode"), "Missing implementation of _encode for %s" % self

        if isinstance(data, (EncodedArray, EncodedRaggedArray)):
            assert data.encoding.is_base_encoding(), "Data is already encoded. " \
                                                     "Can only encode already encoded data if it is base encoded."
            data = data.raw()  # input from here is always "raw"

        if isinstance(data, str):
            out = self._encode_string(data)
        elif isinstance(data, list):
            out = self._encode_list_of_strings(data)
        elif isinstance(data, RaggedArray):
            r = self._ragged_array_as_encoded_array(data)
            assert isinstance(r, (EncodedRaggedArray, RaggedArray))
            return r
        elif isinstance(data, np.ndarray):
            if isinstance(self, NumericEncoding):
                out = self._encode(data)
            else:
                out = EncodedArray(self._encode(data), self)
        else:
            assert False, f"Wrong input type for encode: {type(data)} {data}"

        return out

    def _encode_list_of_strings(self, s: str):
        s = EncodedRaggedArray(
            EncodedArray(
                [ord(c) for ss in s for c in ss], BaseEncoding),
            [len(ss) for ss in s])
        return self._ragged_array_as_encoded_array(s)

    def _ragged_array_as_encoded_array(self, s):
        data = self.encode(s.ravel())
        if isinstance(data, EncodedArray):
            return EncodedRaggedArray(data, s._shape)

        return RaggedArray(data, s._shape)

    def _encode_string(self, string: str):
        s = EncodedArray(np.frombuffer(bytes(string, encoding="ascii"), dtype=np.uint8), BaseEncoding)
        s = self._encode_base_encoded_array(s)
        return s

    def _encode_base_encoded_array(self, encoded_array):
        assert encoded_array.encoding.is_base_encoding()
        encoded_array = self._encode(encoded_array.data)
        if self.is_numeric():
            encoded_array = encoded_array
        else:
            encoded_array = EncodedArray(encoded_array, self)
        return encoded_array

    def decode(self, data):
        if not hasattr(self, "_decode"):
            raise Exception("Missing implementation of _decode for %s" % self)

        if isinstance(data, int):
            return EncodedArray(self._decode(np.atleast_1d(data)), self)
        elif isinstance(data, np.ndarray):
            assert isinstance(self, NumericEncoding), "%s" % data
            return self._decode(data)
        elif isinstance(data, EncodedRaggedArray):
            return EncodedRaggedArray(
                EncodedArray(self._decode(data.raw().ravel()), BaseEncoding), data._shape)
        elif isinstance(data, RaggedArray):
            assert isinstance(self, NumericEncoding), "%s" % data
            return RaggedArray(self._decode(data.ravel()), data._shape)
        elif isinstance(data, EncodedArray):
            return EncodedArray(self._decode(data.raw()), BaseEncoding)
        else:
            raise Exception("Not able to decode %s with %s" % (data, self))

    def is_one_to_one_encoding(self):
        return True


class ASCIIEncoding(OneToOneEncoding):
    def _encode(self, ascii_codes):
        return ascii_codes

    def _decode(self, encoded):
        return encoded

    def __repr__(self):
        return "ASCIIEncoding()"

    def __hash__(self):
        return hash(repr(self))

    def is_base_encoding(self):
        return True

    def __eq__(self, other):
        return isinstance(other, ASCIIEncoding)


class NumericEncoding(OneToOneEncoding):
    def is_numeric(self):
        return True


BaseEncoding = ASCIIEncoding()


def get_base_encodings():
    return [BaseEncoding]  # TODO: add other encodings


class EncodingException(Exception):
    pass


class IncompatibleEncodingsException(Exception):
    pass


class EncodedRaggedArray(RaggedArray):
    """ Class to represnt EnocedArray with different row lengths """

    def __init__(self, data: 'EncodedArray', shape, *args, **kwargs):
        assert isinstance(data, EncodedArray), data
        super().__init__(data.raw(), shape, *args, **kwargs)
        self._encoding = data.encoding

    def copy(self):
        return self.__class__(
            EncodedArray(self.ravel().copy(), self._encoding), self.shape)

    @property
    def _cls(self):
        return lambda data, shape: self.__class__(EncodedArray(data, self._encoding), shape)

    def _set_data_range(self, idx, data):
        super()._set_data_range(idx, as_encoded_array(data, self._encoding).raw())

    def _get_data_range(self, idx):
        return EncodedArray(super()._get_data_range(idx), self._encoding)

    def __repr__(self):
        try:
            return self._proper_repr()
        except Exception as e:
            r = super().__repr__()
            return f'EncodedRaggedArray({r}, {self.encoding()})'

    def _proper_repr(self) -> str:
        if len(self) == 0:
            return ''
        if self.size > 1000:
            rows = [str(row) for row in self[:5]]
        else:
            rows = [f"{row}" for row in self]
        encoding_info = f", {self.encoding}" if not self.encoding.is_base_encoding() else ""
        indent = " " * len("encoded_ragged_array([")
        quotes = "'" if self.encoding.is_one_to_one_encoding() else ""
        lines = [f"{indent}{quotes}{row}{quotes}," for row in rows]
        lines[0] = lines[0].replace(indent, "encoded_ragged_array([")
        if self.size > 1000:
            lines.insert(-1, "...")
        lines[-1] = lines[-1][:-1] + "]" + encoding_info + ")"
        return "\n".join(lines)

    @property
    def encoding(self):
        """Make the encoding of the underlying data avaible"""
        return self._encoding
        # return self.ravel().encoding

    def raw(self):
        return RaggedArray(self.ravel().raw(), self._shape)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """ Convert any data to `EncodedArray` before calling the `ufunc` on them """
        assert isinstance(self.ravel(), EncodedArray), self.ravel()
        inputs = [as_encoded_array(i, self.ravel().encoding).raw() for i in inputs]
        # inputs = _parse_ufunc_inputs(inputs, self.ravel().encoding)
        kwargs = {key: as_encoded_array(val, self.ravel().encoding).raw() for key, val in kwargs.items()}

        ret = super().__array_ufunc__(ufunc, method, *inputs, **kwargs)
        if isinstance(ret.ravel(), EncodedArray):
            return EncodedRaggedArray(ret.ravel(), ret._shape)
        return ret

    def tolist(self):
        return [row.to_string() for row in self]

    def ravel(self):
        return EncodedArray(super().ravel(), self._encoding)


def get_NPSArray(array):
    return array.view(NPSArray)


class EncodedArray(np.lib.mixins.NDArrayOperatorsMixin):
    """ 
    Class for data that could be written as characters, but is represented numpy arrays
    """

    encoding = None

    def __init__(self, data: np.ndarray, encoding):
        """Create an encoded array form raw data and encoding.

        This should seldom be used directly. Only when you for some reason 
        have encoded data that is not yet represented as an `EncodedArray`.

        Use `as_encoded_array` to create encoded arrays from `str` objects



        Parameters
        ----------
        data : np.ndarray
            Raw encoded data
        encoding : Encoding
            The encoding that the data already has

        Examples
        --------
        >>> import bionumpy as bnp
        >>> import numpy as np
        >>> print(EncodedArray(np.array([0, 1, 2, 3]), bnp.DNAEncoding))
        ACGT
        """

        if isinstance(data, EncodedArray):
            assert data.encoding == encoding
            data = data.data
        self.encoding = encoding
        dtype = None if hasattr(data, "dtype") else np.uint8
        # self.data = np.asarray(data, dtype=dtype).view(NPSArray)
        self.data = np.asarray(data, dtype=dtype)
        self.data = get_NPSArray(self.data)
        # assert isinstance(self.data, np.ndarray)

    @property
    def T(self) -> "EncodedArray":
        return self.__class__(self.data.T, self.encoding)

    def copy(self) -> 'EncodedArray':
        return self.__class__(self.data.copy(), self.encoding)

    def __len__(self) -> int:
        return len(self.data)

    def raw(self) -> np.ndarray:
        return self.data.view(np.ndarray)

    def tolist(self) -> str:
        """Converts the data to a string by decoding the data.
        This behaviour is compatible with NumPy's scalar behaviour (only a single element)"""
        return self.to_string()

    def to_string(self) -> str:
        """Converts the data to a string by decoding the data"""
        if not self.encoding.is_one_to_one_encoding():
            return self.encoding.to_string(self.data)
        if hasattr(self.encoding, "_decode"):
            data = self

            raw = self.encoding.decode(data).raw()
            raw = np.atleast_1d(np.asanyarray(raw, dtype=np.uint8))
            text = bytes(raw).decode('ascii')
            return text
            # return "".join([chr(c) for c in self.encoding.decode(data).raw()])
        else:
            data = self.data
            return "".join([chr(c) for c in self.encoding.decode(data)])

    def reshape(self, *args, **kwargs) -> "EncodedArray":
        return self.__class__(self.data.reshape(*args, **kwargs), self.encoding)

    @property
    def size(self) -> int:
        return self.data.size

    @property
    def ndim(self) -> int:
        return self.data.ndim

    @property
    def shape(self) -> Tuple[int]:
        return self.data.shape

    @property
    def strides(self) -> Tuple[int]:
        return self.data.strides

    @property
    def dtype(self) -> np.dtype:
        return self.data.dtype

    def __repr__(self) -> str:
        quotes = "'" if self.encoding.is_one_to_one_encoding() else ""
        if self.encoding.is_base_encoding():
            return f"encoded_array({quotes}{str(self)}{quotes})"
        return f"encoded_array({quotes}{str(self)}{quotes}, {self.encoding})"

    def __str__(self) -> str:
        """Return the data decoded into ASCII string

        Only return the first 20 chars and/or first 20 rows

        Returns
        -------
        str

        """
        if not self.encoding.is_one_to_one_encoding():
            # not possible to decode, get string
            n_dims = len(self.data.shape)
            assert n_dims in [0, 1, 2], "Unsupported number of dimensions for data"

            data_to_show = self.data
            if n_dims == 0:
                return self.encoding.to_string(self.data)
            elif n_dims == 1:
                data_to_show = data_to_show  # [0:20]
            elif n_dims == 2:
                # show first columns and rows
                # raise NotImplemented("Str for n_dims=2 not implemented")
                return self.encoding.to_string(self.data[0:10])  # + "...."
                # data_to_show = None
            text = "[" + ", ".join(self.encoding.to_string(e).strip() for e in data_to_show) + "]"
            return text

        else:
            data = self.data
            if hasattr(self.encoding, "_decode"):
                data = self  # todo: Default after encoding refactoring
                text = self.encoding.decode(data).raw()
            else:
                text = self.encoding.decode(data)

            if len(self.data.shape) == 0:
                return chr(int(text))
            if len(self.data.shape) == 1:
                return "".join(chr(n) for n in text)  # + "..."*(len(text)>20)

            a = np.array([str(self.__class__(seq, self.encoding))
                          for seq in self.data.reshape(-1, self.data.shape[-1])]
                         ).reshape(self.data.shape[:-1])[:20]
            return str(a)

    def __hash__(self) -> int:
        if len(self.shape) <= 1:
            return hash(self.to_string())

    def __getitem__(self, idx) -> "EncodedArray":
        """Delegate the indexing to the underlying numpy array

        Always return EncodedArray object. Even for scalars

        Parameters
        ----------
        idx :
            Any index understood by numpy

        Returns
        -------
        "EncodedArray"

        """
        new_data = self.data.__getitem__(idx)
        if isinstance(new_data, RaggedArray):
            return EncodedRaggedArray(EncodedArray(new_data.ravel(), self.encoding),
                                      new_data._shape)
        return self.__class__(new_data, self.encoding)

    def __setitem__(self, idx, value: "EncodedArray"):
        """ Set the item on the underlying numpy array

        Converts any string values to EncodedArray first

        Parameters
        ----------
        idx : 3
            Any index understood by numpy
        value : "EncodedArray"
            Anything that can be encoded to EncodedArray

        """
        assert isinstance(value, str) or isinstance(value, EncodedArray)

        # if isinstance(value, str):
        #     value = self.encoding.encode(value)
        value = as_encoded_array(value, self.encoding)
        self.data.__setitem__(idx, value.data)

    def __iter__(self):
        return (self.__class__(element, self.encoding) for element in self.data)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Handle numpy ufuncs called on EnocedArray objects

        Only support euqality checks for now. Numeric operations must be performed
        directly on the underlying data.
        """
        if not all(isinstance(a, (
                str, list, EncodedArray,
                EncodedRaggedArray)) for a in inputs):
            return NotImplemented

        if method == "__call__" and ufunc.__name__ in ("equal", "not_equal"):
            inputs = _parse_ufunc_inputs(inputs, self.encoding)
            return ufunc(*inputs)
        return NotImplemented

    def __array_function__(self, func, types, args, kwargs):
        """Handle numpy arrayfunctions called on `EncodedArray` objects

        Limited functionality for now, but can be updated as needed
        """
        if func == np.bincount:
            return np.bincount(args[0].data, *args[1:], **kwargs)
        if func == np.argsort:
            return np.argsort(args[0].data, *args[1:], **kwargs)
        if func == np.concatenate:
            return self.__class__(func([e.data for e in args[0]]), self.encoding)
        if func == np.where:
            return self.__class__(func(args[0], args[1].data, args[2].data), encoding=self.encoding)
        if func == np.zeros_like:
            return self.__class__(func(args[0].data, *args[1:], **kwargs), encoding=self.encoding)
        if func == np.append:
            return self.__class__(func(args[0].data, args[1].data, *args[2:], **kwargs), encoding=self.encoding)
        if func == np.lexsort:
            if not all(issubclass(t, (EncodedArray, np.ndarray)) for t in types):
                return NotImplemented

            args = [a.raw() if isinstance(a, EncodedArray) else np.asarray(a) for a in args[0]]
            return func(args, *kwargs)
        if func == np.full_like:
            return full_like(*args, **kwargs)
        if func == np.insert:
            return self.__class__(func(args[0].data, args[1], args[2].data, *args[3:], **kwargs),
                                  encoding=self.encoding)
        elif func in (np.lib.stride_tricks.sliding_window_view, np.lib.stride_tricks.as_strided):
            return self.__class__(func(args[0].data, *args[1:], **kwargs), self.encoding)

        return NotImplemented
        return super().__array_function__(func, types, args, kwargs)

    def ravel(self) -> "EncodedArray":
        return self.__class__(self.data.ravel(), self.encoding)

    def as_strided(self, *args, **kwargs):
        """Explicitly delegate as_strided, as it is not a arrayfunction

        This should always be used with care, and can lead to some segmentation faults!

        So don't use it :)

        """
        assert isinstance(self.data, np.ndarray) and not np.issubdtype(self.data.dtype, np.object_)
        return self.__class__(np.lib.stride_tricks.as_strided(self.data, *args, **kwargs), self.encoding)


def _parse_ufunc_inputs(inputs, target_encoding):
    for a in inputs:
        assert isinstance(a, (
        str, list, EncodedArray, EncodedRaggedArray)), "%s is not str, list, EncodedArray or EncodedRaggedArray" % type(
            a)
        yield as_encoded_array(a, target_encoding).raw()
        # if isinstance(a, str):
        #     a = target_encoding.encode(a)
        # elif isinstance(a, list):
        #     a = target_encoding.encode(a)
        # else:
        #     # already encoded
        #     pass
        # 
        # yield a.raw()


def _list_of_encoded_arrays_as_encoded_ragged_array(array_list: List[EncodedArray]):
    assert all(isinstance(a, EncodedArray) for a in array_list)
    encoding = array_list[0].encoding
    assert all(a.encoding == encoding for a in array_list)
    if all(a.data.ndim == 0 for a in array_list):
        data = np.array([a.data for a in array_list])
        return EncodedArray(data, encoding)
    data = np.concatenate([a.data for a in array_list])
    shape = [len(a) for a in array_list]
    return EncodedRaggedArray(EncodedArray(data, encoding), shape)


def _is_encoded(data):
    return isinstance(data, (EncodedArray, EncodedRaggedArray))


def bnp_array_func_dispatch(func):
    def new_func(array, *args, **kwargs):
        result = func(array, *args, **kwargs)
        if result is NotImplemented:
            if hasattr(array, '__array_function__'):
                result = array.__array_function__(func, [array] + args, kwargs)
                if result is NotImplemented:
                    raise Exception
        return result


def as_encoded_array(s: object, target_encoding: "Encoding" = None) -> EncodedArray:
    """Main function used to create encoded arrays from e.g. strings orl lists.
    Can be called with already encoded arrays, and will then do nothing.

    If input is `str` or `List[str]` objects, creates an `EncodedArray` or `EncodedRaggedArray`
    object from them with the given encoding.

    If the input is an `EncodedArray` or `EncodedRaggedArray` AND input is BaseEncoded,
    encode the input to the `target_encoding` if possible. If `target_encoding` is None, nothing is done.

    Raw encoded data as `np.ndarray` objects should not be passed to this function. If you have
    already encoded data in `np.ndarray` objects, use the `EncodedArray.__init__`directly

    Parameters
    ----------
    s : str/List[str]/EnocedArray/EncodedRaggedArray
        The data to be represented in an EncodedArray
    target_encoding : Encoding
        The encoding to use in the resulting EncodedArray

    Returns
    -------
    EncodedArray
        Encoded data in an EncodedArray

    default target encoding None:
    if None: encode as base encoding if it is not encoded
    if already encoded: do nothing
    this function is not for changing encoding on stuff

    """
    if isinstance(s, (EncodedArray, EncodedRaggedArray)):
        if target_encoding is None or s.encoding == target_encoding:
            return s
        else:
            if not s.encoding.is_base_encoding():
                if hasattr(s.encoding, 'get_alphabet') and hasattr(target_encoding, 'get_alphabet'):
                    m = s.raw().max()
                    if (s.encoding.get_alphabet()[:m] == target_encoding.get_alphabet()[:m]):
                        if not m < len(target_encoding.get_alphabet()):
                            raise EncodingException(
                                f"Trying to encode already encoded array with encoding {s.encoding} to encoding {target_encoding}. However {s} contains {EncodedArray(m, s.encoding)}")
                        if isinstance(s, EncodedArray):
                            return s.__class__(s.raw(), target_encoding)
                        elif isinstance(s, EncodedRaggedArray):
                            return s.__class__(EncodedArray(s.ravel().raw(), target_encoding), s.shape)
                raise EncodingException("Trying to encode already encoded array with encoding %s to encoding %s. "
                                        "This is not supported. Use the change_encoding function." % (
                                            s.encoding, target_encoding))
    elif target_encoding is None:
        target_encoding = BaseEncoding

    # if numeric encoding and already np-array, this is already encoded
    if target_encoding.is_numeric():
        if type(s) in (np.ndarray, RaggedArray):
            return s
        elif isinstance(s, list) and (len(s)== 0 or isinstance(s[0], (list,Number, np.ndarray))):
            return RaggedArray(s)
    # is already encoded if list and elements are encoded
    elif isinstance(s, list) and len(s) > 0 and isinstance(s[0], EncodedArray):
        return _list_of_encoded_arrays_as_encoded_ragged_array(s)
    if (not isinstance(s, (EncodedArray, EncodedRaggedArray, RaggedArray))) and hasattr(s, 'to_numpy'):
        s = s.to_numpy()
    if isinstance(s, np.ndarray) and (s.dtype == object or np.issubdtype(s.dtype, np.character)):
        s = s.tolist()

    return target_encoding.encode(s)


def full_like(a, fill_value, dtype=None, order='K', subok=True, shape=None):
    assert dtype is None
    assert order == 'K'
    assert subok is True
    fill_value = a.encoding.encode(fill_value)
    return EncodedArray(np.full_like(a.raw(), fill_value, shape=shape), a.encoding)


def from_encoded_array(encoded_array: EncodedArray) -> str:
    """Convert data in an `EncodedArray`/`EncodedRaggedArray into `str`/`List[str]`

    Unlike the `EncodedArray.__str__` this will convert all the data into strings

    Parameters
    ----------
    encoded_array : EncodedArray

    Returns
    -------
    str
        Full string representation

    Examples
    --------
    >>> import bionumpy as bnp
    >>> encoded_array = bnp.DNAEncoding.encode("ACGT")
    >>> print(from_encoded_array(encoded_array))
    ACGT

    >>> encoded_array = bnp.DNAEncoding.encode(["ACGT", "ACGT"])
    >>> print(from_encoded_array(encoded_array))
    ['ACGT', 'ACGT']
    """
    if isinstance(encoded_array, EncodedRaggedArray):
        return [from_encoded_array(row) for row in encoded_array]
    else:
        return "".join(chr(c) for c in encoded_array.encoding.decode(encoded_array).raw())


def change_encoding(encoded_array: Union[EncodedArray, EncodedRaggedArray], new_encoding: Encoding) \
        -> Union[EncodedArray, EncodedRaggedArray]:
    """
    Changes the encoding of an `EncodedArray` or `EncodedRaggedArray` by decoding the data and
    encoding it again with the new encoding.

    Parameters:
    -----------
    encoded_array : EncodedArray/EncodedRaggedArray
        The data to change encoding on
    new_encoding : Encoding
        The new encoding to use

    Returns:
    --------
    EncodedArray/EncodedRaggedArray
        The data with the new encoding

    Examples
    --------
    >>> import bionumpy as bnp
    >>> encoded_array = bnp.as_encoded_array("ACGT", bnp.DNAEncoding)
    >>> encoded_array.raw()
    array([0, 1, 2, 3], dtype=uint8)
    >>> new_encoding = bnp.BaseEncoding
    >>> new_encoded_array = change_encoding(encoded_array, new_encoding)
    >>> new_encoded_array.raw()
    array([65, 67, 71, 84], dtype=uint8)

    """
    assert isinstance(encoded_array, (EncodedArray, EncodedRaggedArray)), \
        "Can only change encoding of EncodedArray or EncodedRaggedArray"

    new_data = new_encoding.encode(
        encoded_array.encoding.decode(encoded_array.ravel())
    )

    if isinstance(encoded_array, EncodedArray):
        return EncodedArray(new_data, new_encoding)
    elif isinstance(encoded_array, EncodedRaggedArray):
        return EncodedRaggedArray(EncodedArray(new_data, new_encoding), encoded_array._shape)


class EncodedLookup:
    def __init__(self, lookup: np.ndarray, encoding: Encoding):
        self._lookup = lookup
        self._encoding = encoding

    def __getitem__(self, key):
        key = self._translate_key(key)
        return self._lookup[key]

    def __setitem__(self, key, value):
        self._lookup[key]=value

    def _translate_key(self, key):
        if isinstance(key, tuple):
            key = tuple(as_encoded_array(i, self._encoding).raw() for i in key)
        else:
            key = as_encoded_array(key, self._encoding).raw()
        return key


def encoded_array_from_nparray(column):
    if hasattr(column, 'raw'):
        column = column.raw()
    if not column.flags['C_CONTIGUOUS']:
        column = column.flatten()
    bytes = column.view(np.uint8).reshape(len(column), -1)
    mask = bytes != 0
    data = bytes[mask]
    return EncodedRaggedArray(EncodedArray(data, BaseEncoding), mask.sum(axis=-1))
