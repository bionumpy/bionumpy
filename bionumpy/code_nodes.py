import numpy as np
from traceback import extract_stack, format_list

class NpDataclassStream:
    def __init__(self, stream, dataclass=None):
        self._stream = stream
        self._dataclass = dataclass
        self.start_buffer = next(self._stream)

    def __getattr__(self, name):
        return LeafNode(self, name)

    def __iter__(self):
        yield self.start_buffer
        yield from self._stream

    def __str__(self):
        return "NpDataclass stream with first buffer:\n" + str(self.start_buffer)


class ComputationNode(np.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, origin_stream, func, args, kwargs=None, stack_trace=None):
        self._origin_stream = origin_stream
        self._func = func
        self._args = args
        self._kwargs = kwargs
        self._buffer = self._do_calc(self._origin_stream.start_buffer)

    def __str__(self):
        return "Array stream with first buffer:\n" + str(self._buffer)

    def ___add__(self, other):
        return ComputationNode(self._origin_stream, np.add, (self, other))

    def __array_ufunc__(self, ufunc, method,  *args, **kwargs):
        stack_trace = "".join(format_list(extract_stack(limit=5))[:-2])
        assert method == "__call__"
        return UfuncNode(self._origin_stream, ufunc, args, kwargs, stack_trace=stack_trace)

    def __array_function__(self, func, types, args, kwargs):
        stack_trace = "".join(format_list(extract_stack(limit=5))[:-2])
        cls = ArrayFuncNode
        if func in reduction_lookup:
            cls = ReductionNode
        return cls(self._origin_stream, func, args, kwargs, stack_trace=stack_trace)

    def _replace_args(self, chunk):
        return [arg._do_calc(chunk) if isinstance(arg, ComputationNode) else arg for arg in self._args]

    def _do_calc(self, chunk):
        return self._func(*[arg._do_calc(chunk) if isinstance(arg, ComputationNode) else arg for arg in self._args])

    def __iter__(self):
        for chunk in self._origin_stream:
            yield self._do_calc(chunk)

def bincount_reduce(bincount_a, bincount_b):
    if bincount_a.size >= bincount_b.size:
        bincount_a[:bincount_b.size] += bincount_b
        return bincount_a
    bincount_b[:bincount_a.size] += bincount_a
    return bincount_b

def mean_reduce(mean_a, mean_b):
    pass

reduction_lookup = {
    np.bincount: bincount_reduce,
    np.mean: mean_reduce}


class UfuncNode(ComputationNode):
    pass


class ArrayFuncNode(ComputationNode):
    pass

class ReductionNode(ComputationNode):
    pass


class LeafNode(ComputationNode):
    def __init__(self, origin_stream, attribute_name):
        self._attribute_name = attribute_name
        self._origin_stream = origin_stream

    def _do_calc(self, chunk):
        return getattr(chunk, self._attribute_name)


def consume(node):
    return list(node)
