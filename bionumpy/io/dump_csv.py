import numpy as np
from numpy.typing import ArrayLike
from typing import List, Dict, Tuple
from .strops import ints_to_strings, int_lists_to_strings, float_to_strings
from ..encoded_array import EncodedArray, Encoding, EncodedRaggedArray, BaseEncoding
from npstructures import RaggedArray
from ..util import is_subclass_or_instance


def dump_csv(data_dict: List[Tuple], sep: str = "\t") -> EncodedArray:
    """Put each field of the dataclass into a column in a buffer.

    Parameters
    data : bnpdataclass
        Data
    """

    funcs = {int: ints_to_strings,
             str: lambda x: x,
             List[int]: int_lists_to_strings,
             float: float_to_strings,
             List[bool]: lambda x: int_lists_to_strings(x.astype(int), sep="")
             }

    def get_func_for_datatype(datatype):
        if is_subclass_or_instance(datatype, Encoding):
            encoding = datatype

            def dynamic(x):
                if isinstance(x, EncodedRaggedArray):
                    # print(repr(x.ravel()))
                    return EncodedRaggedArray(EncodedArray(encoding.decode(x.ravel()), BaseEncoding), x.shape)
                return EncodedArray(encoding.decode(x), BaseEncoding)
            return dynamic
        else:
            return funcs[datatype]

    columns = [get_func_for_datatype(key)(value) for key, value in data_dict]
    # columns = [get_func_for_datatype(field.type)(getattr(data, field.name))
    # for field in dataclasses.fields(data)]

    lengths = np.concatenate([((column.lengths if
                                isinstance(column, RaggedArray)
                                else np.array([
                                            column.shape[-1] if len(column.shape) == 2 else 1
                                              ]*len(column))) + 1
                               )[:, np.newaxis]
                              
                              for column in columns], axis=-1).ravel()
    lines = EncodedRaggedArray(EncodedArray(np.empty(lengths.sum(), dtype=np.uint8), BaseEncoding),
                               lengths)
    n_columns = len(columns)
    for i, column in enumerate(columns):
        if (not isinstance(column, RaggedArray)) and len(column.shape)==1:
            column = EncodedRaggedArray.from_numpy_array(column[:, np.newaxis]) # AND HERE
        lines[i::n_columns, :-1] = column
    lines[:, -1] = sep
    lines[(n_columns - 1)::n_columns, -1] = "\n"
    return lines.ravel()
