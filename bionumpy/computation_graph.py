import numpy as np
from itertools import count, chain
from traceback import extract_stack, format_list
from functools import reduce
import operator


def _add_histograms(histogram_a, histogram_b):
    assert np.all(histogram_a[1] == histogram_b[1]), (histogram_a, histogram_b)
    return ((histogram_a[0]+histogram_b[0]), histogram_a[1])


reductions_map = {
    np.sum: operator.add,
    np.histogram: _add_histograms
    }


class Node(np.lib.mixins.NDArrayOperatorsMixin):
    def _get_buffer(self, i: int):
        return NotImplemented

    def __array_ufunc__(self, ufunc, method,  *args, **kwargs):
        stack_trace = "".join(format_list(extract_stack(limit=5))[:-2])
        assert method == "__call__"
        return ComputationNode(ufunc, args, kwargs, stack_trace=stack_trace)

    def __array_function__(self, func, types, args, kwargs):
        stack_trace = "".join(format_list(extract_stack(limit=5))[:-2])
        comp_node = ComputationNode(func, args, kwargs, stack_trace=stack_trace)
        if func in reductions_map:
            return ReductionNode(comp_node, reductions_map[func])
        return comp_node

    def __str__(self):
        return f'{self.__class__.__name__} with current buffer: {self._current_buffer}'

    def compute(self):
        return NotImplemented

    def get_iter(self):
        for i in count():
            try:
                yield self._get_buffer(i)
            except StopIteration:
                break


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

    def compute(self):
        return np.concataenate(list(self._stream))


class ReductionNode(Node):
    def __init__(self, stream, binary_func):
        self._stream  = stream
        self._binary_func = binary_func

    def compute(self):
        return reduce(self._binary_func, self._stream.get_iter())

    def get_iter(self):
        return accumulate

    @classmethod
    def join(cls, reduction_nodes):
        node = ComputationNode(lambda *args: tuple(args), [node._stream for node in reduction_nodes])
        binary_func = lambda t1, t2: tuple(node._binary_func(e1, e2)
                                           for node, e1, e2 in zip(reduction_nodes, t1, t2))
        return cls(node, binary_func)

    def __str__(self):
        return f'{self._binary_func} reduction of: {self._stream}'


class ComputationNode(Node):
    def __init__(self, func, args, kwargs=None, stack_trace=None):
        self._func = func
        self._args = args
        self._kwargs = kwargs if kwargs is not None else {}
        self._stack_trace = stack_trace
        self._buffer_index = -1
        self._get_buffer(0)

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

    def compute(self):
        return np.concatenate(list(self.get_iter()))


def compute(*args): # func, args, kwargs=None):
    if all(isinstance(a, ReductionNode) for a in args):
        return ReductionNode.join(args).compute()
    else:
        assert not any(isinstance(a, ReductionNode) for a in args)
        node = ComputationNode(lambda *a: tuple(a),
                               args)
        return node.compute()
    # 
    # 
    # func, args = args
    # t = ComputationNode(tuple, args)
    # return t.compute()
    # 
    # 
    # if not any(isinstance(a, Node) for a in args):
    #     return func(*args, **kwargs)
    # if kwargs is None:
    #     kwargs = {}
    # node = ComputationNode(func, args, kwargs)
    # return np.concatenate(list(node.get_iter()))
