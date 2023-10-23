import logging
import types

from .stream import BnpStream

logger = logging.getLogger("streamable")


class streamable:
    """A decorator to make a function applicable to streams of data

    Examples
    --------
    >>> from bionumpy.streams import streamable, BnpStream
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
    def __init__(self, reduction: callable = None):
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
                # logger.debug(f"Running {func.__name__} on chunk {i}")
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
            
            stream_args = [i for i, arg in enumerate(args) if isinstance(arg, (BnpStream, types.GeneratorType))]
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
