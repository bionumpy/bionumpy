import numpy as np


class NpDataclassStream:
    """
    Class to hold a stream/generator of npdataclass objects.
    Works with streamable/bnp_broadcast so that functions decorated
    with those decorator will apply the function to each element in
    a NpDataclassStream. To concatenate all the chunks in the generator
    into one npdataclass, use `NpDataclassStream.join()`
    """
    def __init__(self, stream, buffer_type):
        self._stream = stream
        self._opened = False
        self._peeked = False
        self._buffer = None
        self._buffer_type = buffer_type

    def _peek(self):
        if self._peeked:
            return self._buffer
        self._buffer = next(self._stream, None)
        self._peeked = True
        return self._buffer

    def __iter__(self):
        return self
    
    def __next__(self):
        self._opened = True
        if self._peeked:
            self._peeked = False
            buf = self._buffer
            self._buffer = None
            return buf
        return next(self._stream)

    def __str__(self):
        status = "opened" if self._opened else "unopened"
        return f"""\
{status.capitalize()} stream of {self._buffer_type} buffers. Next buffer:
{self._peek()}
"""

    def element_iter(self):
        return (elem for elem in chunk for chunk in self)

    def join(self):
        return np.concatenate(list(self))


def streamable(func):
    """
    Decorator for applying a function to all chunks of a 
    NpDataclassStream.

    If a NpDataclassStream is provided as an argument to the function,
    the function will be applied to all elements of the stream. 
    Arguments that are not instances of NpDataclassStream will be kept for
    all the fucntion calls
    """

    def _args_stream(args, i):
        args = list(args)
        for arg in args[i]:
            new_args = list(args)
            new_args[i] = arg
            yield new_args

    def new_func(*args, **kwargs):
        stream_args = [i for i, arg in enumerate(args) if isinstance(arg, NpDataclassStream)]
        if len(stream_args) == 0:
            return func(*args, **kwargs)
        
        assert len(stream_args) == 1
        i = stream_args[0]
        args_stream = _args_stream(args, i)
        StreamClass = args[i].__class__
        return StreamClass((func(*new_args, **kwargs) for new_args in args_stream), 
                           buffer_type=args[i]._buffer_type)

    return new_func
