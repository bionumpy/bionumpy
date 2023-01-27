import numpy as np
from itertools import count, chain
from traceback import extract_stack, format_list


class Node(np.lib.mixins.NDArrayOperatorsMixin):
    def _get_buffer(self, i: int):
        return NotImplemented

    def __array_ufunc__(self, ufunc, method,  *args, **kwargs):
        stack_trace = "".join(format_list(extract_stack(limit=5))[:-2])
        assert method == "__call__"
        return ComputationNode(ufunc, args, kwargs, stack_trace=stack_trace)

    def __array_function__(self, func, types, args, kwargs):
        stack_trace = "".join(format_list(extract_stack(limit=5))[:-2])
        return ComputationNode(func, args, kwargs, stack_trace=stack_trace)

    def __str__(self):
        return f'{self.__class__.__name__} with current buffer: {self._current_buffer}'


class StreamNode(Node):
    def __init__(self, stream):
        self._stream = stream
        self._current_buffer = None
        self._buffer_index = -1
        self._get_buffer(0)

    def _get_buffer(self, i: int):
        assert self._buffer_index in (i, i-1),  (i, self._buffer_index)
        if i > self._buffer_index:
            self._current_buffer = next(self._stream)
            self._buffer_index += 1
        return self._current_buffer

class ComputationNode(Node):
    def __init__(self, func, args, kwargs=None, stack_trace=None):
        self._func = func
        self._args = args
        self._kwargs = kwargs if kwargs is not None else {}
        self._stack_trace = stack_trace
        self._buffer_index = -1
        self._get_buffer(0)

    def __len__(self):
        return NotImplemented

    def __getitem__(self, item):
        return ComputationNode(self.__class__,
                               (self, item))

    def max(self, axis=None, **kwargs):
        assert axis == -1, axis
        return np.max(self, axis=-1, **kwargs)

    def mean(self, axis=None):
        assert axis == -1, axis
        return np.mean(self, axis=-1)

    def _get_buffer(self, i: int):
        assert self._buffer_index in (i, i-1),  (i, self._buffer_index)
        if i > self._buffer_index:
            args = [a._get_buffer(i) if isinstance(a, Node)
                    else a for a in self._args]
            kwargs = {key: (v._get_buffer(i) if isinstance(v, Node) else v)
                      for key, v in self._kwargs.items()}
            self._current_buffer = self._func(*args, **kwargs)
            self._buffer_index += 1
        return self._current_buffer

    def __iter__(self):
        for i in count():
            try:
                yield self._get_buffer(i)
            except StopIteration:
                break

def compute(func, args, kwargs=None):
    if kwargs is None:
        kwargs = {}
    node = ComputationNode(func, args, kwargs)
    return np.concatenate(list(node))
    
    
