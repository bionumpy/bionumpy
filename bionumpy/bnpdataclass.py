import dataclasses
from typing import List
from npstructures.npdataclasses import npdataclass
from npstructures import RaggedArray
import numpy as np
from .encoded_array import EncodedArray, as_encoded_array, EncodedRaggedArray
from .encodings import Encoding, NumericEncoding
from .util import is_subclass_or_instance

def bnpdataclass(base_class):
    def _implicit_format_conversion(cls, obj):
        for field in dataclasses.fields(obj):
            pre_val = getattr(obj, field.name)
            if field.type in (int, float, -1):
                val = np.asanyarray(pre_val)
            elif field.type == str:
                assert isinstance(pre_val, (str, list, EncodedArray, EncodedRaggedArray, RaggedArray)), (field, pre_val)
                val = as_encoded_array(pre_val)
            elif is_subclass_or_instance(field.type, Encoding):
                if is_subclass_or_instance(field.type, NumericEncoding):
                    assert isinstance(pre_val, (str, list, EncodedArray, EncodedRaggedArray, RaggedArray, np.ndarray)), (field, pre_val)
                else:
                    assert isinstance(pre_val, (str, list, EncodedArray, EncodedRaggedArray)), (field, pre_val)
                val = as_encoded_array(pre_val, field.type)
            elif field.type == List[int]:
                val = RaggedArray(pre_val) if not isinstance(pre_val, RaggedArray) else pre_val
            else:
                assert False, field.type
            setattr(obj, field.name, val)

    new_class = npdataclass(base_class)
    new_class._implicit_format_conversion = classmethod(_implicit_format_conversion)
    return new_class

