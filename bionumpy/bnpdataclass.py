import dataclasses
from typing import List
from npstructures.npdataclasses import npdataclass
from npstructures import RaggedArray
import numpy as np
from .sequences import as_sequence_array, as_encoded_sequence_array, EncodedArray
from .encodings import Encoding


def bnpdataclass(base_class):
    def _implicit_format_conversion(cls, obj):
        for field in dataclasses.fields(obj):
            pre_val = getattr(obj, field.name)
            if field.type in (int, float, -1):
                val = np.asanyarray(pre_val)
            elif field.type == str:
                val = as_sequence_array(pre_val)
            elif (isinstance(field.type, type) and issubclass(field.type, (EncodedArray, Encoding))) or isinstance(field.type, Encoding):
                val = as_encoded_sequence_array(pre_val, field.type)
            elif field.type == List[int]:
                val = RaggedArray(pre_val) if not isinstance(pre_val, RaggedArray) else pre_val
            else:
                assert False, field.type
            setattr(obj, field.name, val)

    new_class = npdataclass(base_class)
    new_class._implicit_format_conversion = classmethod(_implicit_format_conversion)
    return new_class

