import numpy as np
import logging
from npstructures import RaggedArray
import dataclasses

logger = logging.getLogger(__name__)


def as_strided(arr, *args, **kwargs):
    if hasattr(arr, "as_strided"):
        return arr.as_strided(*args, **kwargs)
    assert not np.issubdtype(arr.dtype, np.object_), arr
    return np.lib.stride_tricks.as_strided(arr, *args, **kwargs)


def convolution(func):
    def new_func(_sequence, window_size, *args, **kwargs):
        shape, sequence = (_sequence.shape, _sequence.ravel())
        convoluted = func(sequence, window_size, *args, **kwargs)
        if isinstance(shape, tuple):
            out = as_strided(convoluted, shape)
        else:
            out = RaggedArray(convoluted, shape)

        return out[..., : (-window_size + 1)]

    return new_func


def rolling_window_function(func):
    def new_func(_sequence, window_size, *args, **kwargs):
        shape, sequence = (_sequence.shape, _sequence.ravel())
        windows = np.lib.stride_tricks.sliding_window_view(sequence, window_size)
        convoluted = func(windows, window_size, *args, **kwargs)
        if isinstance(_sequence, RaggedArray):
            out = RaggedArray(convoluted, shape)
        elif isinstance(_sequence, np.ndarray):
            out = as_strided(convoluted, shape)
        return out[..., : (-window_size + 1)]

    return new_func


def pprint_one(sequence):
    return "".join(chr(c) for c in sequence)


def pprint(sequences):
    if isinstance(sequences, RaggedArray):
        return [pprint_one(s) for s in sequences]
    elif isinstance(sequences, np.ndarray):
        if len(sequences.shape) == 1:
            return pprint_one(sequences)
        return [pprint(s) for s in sequences]


def plot(obj):
    if not hasattr(obj, "__plot__"):
        logger.warning(f"{obj} has no __plot__ method")


def apply_to_npdataclass(attribute_name):
    def decorator(func):
        def new_func(np_dataclass, *args, **kwargs):
            if not isinstance(np_dataclass, NpDataclass):
                return func(np_dataclass)
            
            return dataclasses.replace(**{attribute_name: func(getattr(np_dataclass, attribute_name), *args, **kwargs)})
        return new_func
    return decorator

def is_subclass_or_instance(obj, c):
    return (isinstance(obj, type) and issubclass(obj, c)) or isinstance(obj, c)

