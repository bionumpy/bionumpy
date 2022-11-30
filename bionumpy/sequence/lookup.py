from ..encodings import BaseEncoding
from ..encoded_array import as_encoded_array


class Lookup:
    def __init__(self, values, encoding=BaseEncoding):
        self._values = values
        self._encoding = encoding

    def __getitem__(self, idx):
        if isinstance(idx, tuple):
            idx = tuple(as_encoded_array(i, self._encoding).raw() if not (isinstance(i, slice) or i is Ellipsis) else i
                        for i in idx)
        else:
            idx = as_encoded_array(idx, self._encoding).raw() if not (isinstance(idx, slice) or idx is Ellipsis) else idx
        return self._values[idx]

    def __setitem__(self, idx, value):
        if isinstance(idx, tuple):
            idx = tuple(as_encoded_array(i, self._encoding).raw() for i in idx)
        else:
            idx = as_encoded_array(idx, self._encoding).raw()
        return self._values.__setitem__(idx, value)
