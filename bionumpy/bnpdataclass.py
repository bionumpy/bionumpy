import dataclasses
from npstructures.npdataclasses import npdataclass
import numpy as np
from .sequences import as_sequence_array
from .encodings import Encoding


def bnpdataclass(base_class):
    def _implicit_format_conversion(cls, obj):
        for field in dataclasses.fields(obj):
            pre_val = getattr(obj, field.name)
            if field.type in (int, float):
                val = np.asanyarray(pre_val)
            elif field.type == str:
                val = as_sequence_array(pre_val)
            elif issubclass(field.type, Encoding):
                val = as_encoded_sequence_array(pre_val, field.type)
            setattr(obj, field.name, val)

    def _implicit_format_conversion_single(self):
        print("IS_CALLED")
        for field in dataclasses.fields(self):
            pre_val = getattr(self, field.name)
            if field.type == str:
                val = as_sequence_array(pre_val)
            elif issubclass(field.type, Encoding):
                val = as_encoded_sequence_array(pre_val, field.type)
            else:
                val = pre_val
            setattr(self, field.name, val)
    

    new_class = npdataclass(base_class)
    new_class._implicit_format_conversion = classmethod(_implicit_format_conversion)
    return new_class
