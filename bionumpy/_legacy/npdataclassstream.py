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

class streamable:
    """A decorator to make a function applicable to streams of data

    Examples
    --------

    >>> @streamable(list)
    ... def add(a, b):
    ...       return a + b

    >>> add(10, 15)
    25
    >>> stream = BnpStream(iter(range(10)))
    >>> add(stream, 15)
    [15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

    >>> @streamable(sum)
    ... def mul(a, b):
    ...       return a * b
    >>> stream = BnpStream(iter([3, 4]))
    >>> mul(stream, 2)
    14
"""

    def __init__(self, reduction: callable=None):
        """Take an optional reduction operator that will be called on the result stream

        Parameters
        ----------
        reduction : callable
            A reduction function that can be called on a stream of
            data
        """
        
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

    def __call__(self, func: callable) -> callable:
        """Return a new function that applies the input function to all chunks in a stream

        Parameters
        ----------
        func : callable
            A normal function that operates on one chunk of data

        Returns
        -------
        callable
            A new function that applies the original function to every
            chunk in a stream of data

        Examples
        --------
        FIXME: Add docs.

        """
        
        def log(sequence):
            for i, args in enumerate(sequence):
                logger.info(f"Running {func.__name__} on chunk {i}")
                yield args

        def new_func(*args, **kwargs):
            """Apply a function to all chunks in a stream if any of the arguments are `BnpStream`

            If no arguments are BnpStream, then the function is simply
            called on the given `args` and `kwargs`. If however one of
            the arguments is a `BnpStream` the the function is applied
            to all chunks in the stream. If the `reduction` parameter
            was given, then the reduction is also called on the
            resulting stream.

            Parameters
            ----------
            *args :
            **kwargs :

            """
            
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


def quantile(array, quantiles, axis=None):
    """
    """
    hist = bincount(array)
    cumulative = np.cumsum(hist)
    idxs = np.searchsorted(cumulative, quantiles*cumulative[-1])
    return idxs
