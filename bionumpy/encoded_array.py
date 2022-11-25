from npstructures import RaggedArray
from npstructures.mixin import NPSArray
from typing import Tuple
import numpy as np


class EncodingException(Exception):
    pass


class IncompatibleEncodingsException(Exception):
    pass


class EncodedRaggedArray(RaggedArray):
    """ Class to represnt EnocedArray with different row lengths """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert isinstance(self._data, EncodedArray)

    def __repr__(self) -> str:
        if self.size>1000:
            rows = [str(row) for row in self[:5]]
        else:
            rows = [f"{row}" for row in self]
        encoding_info = f", {self.encoding}" if not self.encoding.is_base_encoding()  else ""
        indent = " "*len("encoded_ragged_array([")
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
        return self._data.encoding

    def raw(self):
        return RaggedArray(self._data.raw(), self.shape)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        from .encoded_array_functions import as_encoded_array
        """ Convert any data to `EncodedArray` before calling the `ufunc` on them """
        assert isinstance(self._data, EncodedArray), self._data
        inputs = [as_encoded_array(i, self._data.encoding).raw() for i in inputs]
        #inputs = _parse_ufunc_inputs(inputs, self._data.encoding)
        kwargs = {key: as_encoded_array(val, self._data.encoding).raw() for key, val in kwargs.items()}

        ret = super().__array_ufunc__(ufunc, method, *inputs, **kwargs)
        if isinstance(ret._data, EncodedArray):
            return EncodedRaggedArray(ret._data, ret.shape)
        return ret

    def tolist(self):
        return [row.to_string() for row in self]


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

        """
        if isinstance(data, EncodedArray):
            assert data.encoding == encoding
            data = data.data
        self.encoding = encoding
        dtype = None if hasattr(data, "dtype") else np.uint8
        #self.data = np.asarray(data, dtype=dtype).view(NPSArray)
        self.data = np.asarray(data, dtype=dtype)
        self.data = get_NPSArray(self.data)
        #assert isinstance(self.data, np.ndarray)

    def __len__(self) -> int:
        return len(self.data)

    def raw(self) -> np.ndarray:
        return self.data.view(np.ndarray)

    def to_string(self) -> str:
        if hasattr(self, "_decode"):
            # new system, can be used in all cases after refactoring
            data = self
        else:
            data = self.data
        return "".join([chr(c) for c in self.encoding.decode(data)])

    def reshape(self, *args, **kwargs) -> "EncodedArray":
        return self.__class__(self.data.reshape(*args, **kwargs), self.encoding)

    @property
    def size(self) -> int:
        return self.data.size

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
                data_to_show = data_to_show[0:20]
            elif n_dims == 2:
                # show first columns and rows
                #raise NotImplemented("Str for n_dims=2 not implemented")
                return self.encoding.to_string(self.data)[0:40] + "..."
                # data_to_show = None
            text = "[" + ", ".join(self.encoding.to_string(e) for e in data_to_show) + "]"
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
                return "".join(chr(n) for n in text[:20]) + "..."*(len(text)>20)

            a = np.array([str(self.__class__(seq, self.encoding))
                          for seq in self.data.reshape(-1, self.data.shape[-1])]
                         ).reshape(self.data.shape[:-1])[:20]
            return str(a)


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
                                      new_data.shape)
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

        if isinstance(value, str):
            from .encoded_array_functions import encode_string
            value = encode_string(value, self.encoding)

        self.data.__setitem__(idx, value.data)

    def __iter__(self):
        return (self.__class__(element, self.encoding) for element in self.data)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Handle numpy ufuncs called on EnocedArray objects

        Only support euqality checks for now. Numeric operations must be performed
        directly on the underlying data.
        """

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
        if func == np.concatenate:
            return self.__class__(func([e.data for e in args[0]]), self.encoding)
        if func == np.where:
            return self.__class__(func(args[0], args[1].data, args[2].data), encoding = self.encoding)
        if func == np.zeros_like:
            return self.__class__(func(args[0].data, *args[1:], **kwargs), encoding=self.encoding)
        if func == np.append:
            return self.__class__(func(args[0].data, args[1].data, *args[2:], **kwargs), encoding=self.encoding)
        if func == np.insert:
            return self.__class__(func(args[0].data, args[1], args[2].data, *args[3:], **kwargs), encoding = self.encoding)
        elif func in (np.lib.stride_tricks.sliding_window_view, np.lib.stride_tricks.as_strided):
            return self.__class__(func(args[0].data, *args[1:], **kwargs), self.encoding)
        
        return NotImplemented
        return super().__array_function__(func, types, args, kwargs)

    def ravel(self):
        return self.__class__(self.data.ravel(), self.encoding)

    def as_strided(self, *args, **kwargs):
        """Explicitly delegate as_strided, as it is not a arrayfunction

        This should always be used with care, and can lead to some segmentation faults!

        So don't use it :)

        """
        assert isinstance(self.data, np.ndarray) and not np.issubdtype(self.data.dtype, np.object_)
        return self.__class__(np.lib.stride_tricks.as_strided(self.data, *args, **kwargs), self.encoding)


def _parse_ufunc_inputs(inputs, target_encoding):
    from .encoded_array_functions import encode_string, encode_list_of_strings
    for a in inputs:
        assert isinstance(a, (str, list, EncodedArray, EncodedRaggedArray)), repr(a)
        if isinstance(a, str):
            a = encode_string(a, target_encoding)
        elif isinstance(a, list):
            a = encode_list_of_strings(a, target_encoding)
        yield a.raw()



