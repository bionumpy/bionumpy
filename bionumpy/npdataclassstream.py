import numpy as np
from functools import reduce
import dataclasses
import logging
logger = logging.getLogger("streamable")


class BnpStream:
    def __init__(self, stream, first_buffer=None):
        self._stream = stream
        self._first_buffer = first_buffer

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._stream)


class ArrayStream(BnpStream):
    def __init__(self, stream, first_buffer=None):
        self._stream = stream
        self._first_buffer = first_buffer

    def ravel(self):
        return ArrayStream(arr.ravel() for arr in self)


class NpDataclassStream(BnpStream):
    """
    Class to hold a stream/generator of npdataclass objects.
    Works with streamable/bnp_broadcast so that functions decorated
    with those decorator will apply the function to each element in
    a NpDataclassStream. To concatenate all the chunks in the generator
    into one npdataclass, use `NpDataclassStream.join()`
    """
    def __init__(self, stream, dataclass):
        self._stream = stream
        self._opened = False
        self._peeked = False
        self._buffer = None
        self._dataclass = dataclass

    def _peek(self):
        if self._peeked:
            return self._buffer
        self._buffer = next(self._stream, None)
        self._peeked = True
        return self._buffer
        return next(self._stream)

    def __iter__(self):
        return self
    
    def __next__(self):
        return next(self._stream)

        self._opened = True
        if self._peeked:
            self._peeked = False
            buf = self._buffer
            self._buffer = None
            return buf

    def __str__(self):
        status = "opened" if self._opened else "unopened"
        return f"""\
{status.capitalize()} stream of {self._buffer_type} buffers. Next buffer:
{self._peek()}
"""

    def __getattr__(self, attribute_name):
        if not attribute_name in {f.name for f in dataclasses.fields(self._dataclass)}:
            raise Exception(f"{self._dataclass} has no attribute {attribute_name}")
        return ArrayStream(getattr(chunk, attribute_name) for chunk in self)
                            

    def element_iter(self):
        return (elem for elem in chunk for chunk in self)

    def join(self):
        return np.concatenate(list(self))


class streamable:
    def __init__(self, reduction=None):
        self._reduction = reduction

    @staticmethod
    def _args_stream(args, stream_indices):
        args = list(args)
        streams = tuple(args[i] for i in stream_indices)
        for stream_args in zip(*streams):
            new_args = list(args)
            for i, stream_arg in zip(stream_indices, stream_args):
                new_args[i] = stream_arg
            yield new_args

    def __call__(self, func):

        def log(sequence):
            for i, args in enumerate(sequence):
                logger.info(f"Running {func.__name__} on chunk {i}")
                yield args

        def new_func(*args, **kwargs):

            stream_args = [i for i, arg in enumerate(args) if isinstance(arg, BnpStream)]
            if len(stream_args) == 0:
                return func(*args, **kwargs)

            args_stream = log(self._args_stream(args, stream_args))
            
            stream = (func(*new_args, **kwargs) for new_args in args_stream)
            if self._reduction is None:
                return BnpStream(stream)
            else:
                return self._reduction(stream)
            StreamClass = args[stream_args[0]].__class__
            return StreamClass((func(*new_args, **kwargs) for new_args in args_stream),
                               dataclass=args[stream_args[0]].dataclass)

        return new_func


def bincount_reduce(bincount_a, bincount_b):
    if bincount_a.size >= bincount_b.size:
        bincount_a[:bincount_b.size] += bincount_b
        return bincount_a
    bincount_b[:bincount_a.size] += bincount_a
    return bincount_b


bincount = streamable(lambda x: reduce(bincount_reduce, x))(np.bincount)


def histogram_reduce(histograms):
    hist, edge = next(histograms)
    hist = sum((h[0] for h in histograms))+hist
    return hist, edge


histogram = streamable(histogram_reduce)(np.histogram)


@streamable(sum)
def sum_and_n(array, axis=None):
    if axis is None:
        n = array.size
    elif axis == 0:
        n = len(array)
    return np.append(np.sum(array, axis=axis), n)


@streamable()
def _rowmean(array, axis=None):
    return np.mean(array, axis=axis)


def mean(array, axis=None):
    if (axis is not None) and axis != 0:
        return _rowmean(array, axis)
    t = sum_and_n(array, axis=axis)
    return t[:-1]/t[-1]
