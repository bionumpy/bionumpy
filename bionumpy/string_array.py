import numpy as np

from .encoded_array import EncodedRaggedArray, BaseEncoding


class StringArray(np.lib.mixins.NDArrayOperatorsMixin):
    wrapped_functions = ['size', 'shape', 'ndim', '__len__']

    def __init__(self, data):
        self._data = np.asanyarray(data, dtype='S')

    def __repr__(self):
        if self._data.ndim>=1:
            return repr(self._data[:10])
        else:
            return self._data.__repr__()

    def raw(self):
        return self._data

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Handle numpy ufuncs called on EnocedArray objects

        Only support euqality checks for now. Numeric operations must be performed
        directly on the underlying data.
        """

        if method == "__call__" and ufunc.__name__ in ("equal", "not_equal"):
            inputs = [self._convert_input(input) for input in inputs]
            assert len(inputs) == 2, 'only binary operations are supported'
            if ufunc.__name__ == 'equal':
                return inputs[0] == inputs[1]
            elif ufunc.__name__ == 'not_equal':
                return inputs[0] != inputs[1]
            return ufunc(*inputs)
        return NotImplemented

    def __getitem__(self, item):
        return self.__class__(self._data[item])

    def __setitem__(self, item, value):
        self._data[item] == self._convert_input(item, value)

    def __getattr__(self, name):
        if name in self.wrapped_functions:
            return getattr(self._data, name)
        raise AttributeError(f'{self.__class__.__name__} has no attribute {name}')

    def _convert_input(self, value):
        if isinstance(value, str):
            return bytes(value)
        elif isinstance(value, self.__class__):
            return value.raw()
        return np.asanyarray(value, dtype='S')

    def __len__(self):
        return len(self._data)

def printio(func):
    def new_func(*args, **kwargs):
        print(args)
        res = func(*args, **kwargs)
        print('>', res)
        return res
    return new_func

# @printio
def string_array(input_data):
    if isinstance(input_data, list) and len(input_data) > 0 and isinstance(input_data[0], StringArray):
        return string_array([i.raw() for i in input_data])
    if isinstance(input_data, (list, str)):
        return StringArray(np.array(input_data, dtype='S'))
    elif isinstance(input_data, np.ndarray):
        return StringArray(input_data)
    elif isinstance(input_data, StringArray):
        return input_data.copy()
    elif isinstance(input_data, EncodedRaggedArray):
        assert input_data.encoding == BaseEncoding, input_data.encoding
        array = input_data.raw().as_padded_matrix(side='right')
        n_bytes = array.shape[-1]
        return StringArray(array.ravel().view(f'|S{n_bytes}'))
    else:
        raise TypeError(f'Cannot convert {input_data} to StringArray ({type(input_data)})')


def as_string_array(input_data):
    if isinstance(input_data, StringArray):
        return input_data
    return string_array(input_data)