import numpy as np
from itertools import count, chain
from traceback import extract_stack, format_list
from functools import reduce
from abc import ABC, abstractmethod
import operator


class ComputationException(Exception):
    pass


def _add_histograms(histogram_a, histogram_b):
    assert np.all(histogram_a[1] == histogram_b[1]), (histogram_a, histogram_b)
    return ((histogram_a[0]+histogram_b[0]), histogram_a[1])


def mean_reduction(sum_and_n_a, sum_and_n_b):
    return (sum_and_n_a[0]+sum_and_n_b[0],
            sum_and_n_a[1]+sum_and_n_b[1])


def sum_and_n(array, axis=None):
    if array.size == 0:
        return 0, 0

    s = np.sum(array, axis=axis)
    assert axis in (None, 0, -2), axis
    if axis is None:
        n = array.size
    elif axis == 0:
        if hasattr(array.shape[-1], '__len__'):
            if hasattr(array, 'col_counts'):
                n = array.col_counts()
            else:
                n = np.bincount(array.shape[-1])
                n = np.cumsum(n[::-1])[::-1][1:]
        else:
            n = len(array)
    return (s, n)

reductions_map = {
    np.sum: operator.add,
    np.histogram: _add_histograms,
    }


class Node(np.lib.mixins.NDArrayOperatorsMixin, ABC):

    @abstractmethod
    def _get_buffer(self, i: int):
        return NotImplemented

    def __array_ufunc__(self, ufunc, method,  *args, **kwargs):
        stack_trace = "".join(format_list(extract_stack(limit=5)))
        assert method == "__call__"
        return ComputationNode(ufunc, args, kwargs, stack_trace=stack_trace)

    def __array_function__(self, func, types, args, kwargs):
        stack_trace = "".join(format_list(extract_stack(limit=10))[:-2])
        if func == np.mean and ('axis' not in kwargs or kwargs['axis'] in (None, 0)):
            comp_node = ComputationNode(sum_and_n, args, kwargs, stack_trace=stack_trace)
            return ReductionNode(comp_node, mean_reduction, lambda sn: sn[0]/sn[1])
        else:
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
        return np.concatenate(list(self._stream))


class ReductionNode(Node):
    def __init__(self, stream, binary_func, post_process=None):
        self._stream  = stream
        self._binary_func = binary_func
        self._post_process = post_process

    def _get_buffer(self, i: int):
        raise NotImplementedError

    def compute(self):
        r = reduce(self._binary_func, self._stream.get_iter())
        if self._post_process is not None:
            r = self._post_process(r)
        return r

    def get_iter(self):
        return accumulate

    @classmethod
    def join(cls, reduction_nodes):
        node = ComputationNode(lambda *args: tuple(args), [node._stream for node in reduction_nodes])
        binary_func = lambda t1, t2: tuple(node._binary_func(e1, e2)
                                           for node, e1, e2 in zip(reduction_nodes, t1, t2))
        post_process = lambda t: (e if node._post_process is None else node._post_process(e) for e, node in zip(t, reduction_nodes))
        return cls(node, binary_func, post_process)

    def __str__(self):
        return f'{self._binary_func} reduction of: {self._stream}'


class ComputationNode(Node):
    def __init__(self, func, args, kwargs=None, stack_trace=None):
        self._func = func
        self._args = args
        self._kwargs = kwargs if kwargs is not None else {}
        self._stack_trace = stack_trace
        self._buffer_index = -1
        self._stack_trace = "".join(format_list(extract_stack(limit=5))[:-2])
        self._get_buffer(0)

    def __getitem__(self, item):
        return ComputationNode(lambda obj, item: obj[item],
                               (self, item))

    def max(self, axis=None, **kwargs):
        assert axis == -1, axis
        return np.max(self, axis=-1, **kwargs)

    def mean(self, axis=None):
        if axis  in (-1, 1, 0):
            return np.mean(self, axis=axis)
        raise ValueError('invalid axis for mean', axis)

    def sum(self, *args, **kwargs):
        return np.sum(self, *args, **kwargs)

    def _get_buffer(self, i: int):
        assert self._buffer_index in (i, i-1),  (i, self._buffer_index)
        if i <= self._buffer_index:
            return self._current_buffer
        args = [a._get_buffer(i) if isinstance(a, Node)
                else a for a in self._args]
        kwargs = {key: (v._get_buffer(i) if isinstance(v, Node) else v)
                  for key, v in self._kwargs.items()}
        try:
            self._current_buffer = self._func(*args, **kwargs)
        except Exception as e:
            if isinstance(e, ComputationException):
                raise
            raise ComputationException(f'Error in computation of:\n {self._stack_trace}') from e

        self._buffer_index += 1
        
        return self._current_buffer

    def compute(self):
        return np.concatenate(list(self.get_iter()))


class JoinNode(ComputationNode):
    
    def compute(self):
        buffer_list = None
        for buffer_tuple in self.get_iter():
            if buffer_list is None:
                buffer_list = [list() for _ in buffer_tuple]
            for l, b in zip(buffer_list, buffer_tuple):
                l.append(b)
        return [np.concatenate(column) for column in buffer_list]
    # return [np.concatenate(column) for column in zip(*self.get_iter())]


def _compute(*args):
    if not any(isinstance(a, Node) for a in args):
        return args
    if all(isinstance(a, ReductionNode) for a in args):
        return ReductionNode.join(args).compute()
    else:
        assert not any(isinstance(a, ReductionNode) for a in args)
        node_idxs = [i for i, a in enumerate(args) if isinstance(a, Node)]
        results = JoinNode(lambda *a: tuple(a),
                           [args[i] for i in node_idxs]).compute()
        args = list(args)
        for i, idx in enumerate(node_idxs):
            args[idx] = results[i]
        return args


def compute(args):
    if isinstance(args, dict):
        return dict(zip(args.keys(), _compute(*args.values())))
    elif isinstance(args, (list, tuple)):
        return _compute(*args)
    elif isinstance(args, Node):
        return args.compute()
    return args
    

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
