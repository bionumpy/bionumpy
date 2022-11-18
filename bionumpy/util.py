import numpy as np
import logging
from npstructures import RaggedArray
import dataclasses
import sys

logger = logging.getLogger(__name__)


def as_strided(arr, shape=None, strides=None, **kwargs):
    logger.warning(f"{(arr, shape, strides)}")
    sys.stdout.flush()
    sys.stderr.flush()
    # assert strides is not None, (arr, shape, strides)
    if strides is None:
        assert len(arr.shape) == 1
        if len(shape) == 2:
            strides = (shape[-1]*arr.strides[-1], arr.strides[-1])
        elif len(shape) == 1:
            strides = (arr.strides[-1],)
        else:
            assert False, (arr, shape, len(shape))

    if hasattr(arr, "as_strided"):
        return arr.as_strided(shape, strides, **kwargs)
    assert not np.issubdtype(arr.dtype, np.object_), arr
    return np.lib.stride_tricks.as_strided(arr, shape, strides, **kwargs)


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

