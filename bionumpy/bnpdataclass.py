import dataclasses
from npstructures.npdataclasses import npdataclass
import numpy as np
from .sequences import as_sequence_array
from .encodings import Encoding


def bnpdataclass(base_class):
    def _implicit_format_conversion(self):
        for field in dataclasses.fields(self):
            pre_val = getattr(self, field.name)
            if field.type in (int, float):
                val = np.asanyarray(pre_val)
            elif field.type == str:
                val = as_sequence_array(pre_val)
            elif issubclass(field.type, Encoding):
                val = as_encoded_sequence_array(pre_val, field.type)
            setattr(self, field.name, val)
    new_class = npdataclass(base_class)
    new_class._implicit_format_conversion = _implicit_format_conversion
    return new_class
