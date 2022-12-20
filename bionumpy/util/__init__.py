import numpy as np
import numpy.typing as npt
import logging
from npstructures import RaggedArray
import dataclasses
import sys
from . import testing
from . import typing
logger = logging.getLogger(__name__)


def as_strided(arr, shape=None, strides=None, **kwargs):
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


def is_subclass_or_instance(obj, c):
    return (isinstance(obj, type) and issubclass(obj, c)) or isinstance(obj, c)

def interleave(array_a: npt.ArrayLike, array_b: npt.ArrayLike) -> npt.ArrayLike:
    c = np.empty((array_a.size + array_b.size,), dtype=array_a.dtype)
    c[0::2] = array_a
    c[1::2] = array_b
    return c
