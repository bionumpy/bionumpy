import numpy as np
from npstructures import RaggedArray

from ..encoded_array import EncodedArray, EncodedRaggedArray


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
            return self.pd.Series(object.tolist(), dtype='string')
        if isinstance(object, RaggedArray):
            return list(object)
        return object.todict()


pandas_adaptor = PandasAdaptor()