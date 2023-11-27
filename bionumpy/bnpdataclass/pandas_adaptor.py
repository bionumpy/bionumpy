import numpy as np
from npstructures import RaggedArray

from ..encoded_array import EncodedArray, EncodedRaggedArray
from ..string_array import StringArray


class PandasAdaptor:
    def __init__(self):
        self._pd = None

    @property
    def pd(self):
        if self._pd is None:
            import pandas as pd
            self._pd = pd
        return self._pd

    def get_data_frame(self, dict_object: dict):
        return self.pd.DataFrame(dict_object)

    def pandas_converter(self, object):
        if isinstance(object, np.ndarray):
            return object
        if isinstance(object, (EncodedArray, EncodedRaggedArray)):
            l = object.tolist()
            if isinstance(l, str):
                l = list(l)
            return self.pd.Series(l, dtype='string')
        if isinstance(object, RaggedArray):
            return list(object)
        if isinstance(object, StringArray):
            return object.tolist()
        return object.todict()


pandas_adaptor = PandasAdaptor()