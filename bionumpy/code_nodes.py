import numpy as np


class Mockup:
    def __init__(self, stream):
        self._stream = stream

    def __getattr__(self, name):
        return LeafNode(self._stream, name)


class ComputationNode(np.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, origin_stream, func, args, kwargs=None):
        self._origin_stream = origin_stream
        self._func = func
        self._args = args
        self._kwargs = kwargs

    def ___add__(self, other):
        return ComputationNode(self._origin_stream, np.add, (self, other))

    def __array_ufunc__(self, ufunc, method,  *args, **kwargs):
        assert method == "__call__"
        return UfuncNode(self._origin_stream, ufunc, args, kwargs)

    def __array_function__(self, func, types, args, kwargs):
        return ArrayFuncNode(self._origin_stream, func, args, kwargs)

    def _replace_args(self, chunk):
        return [arg._do_calc(chunk) if isinstance(arg, ComputationNode) else arg for arg in self._args]

    def _do_calc(self, chunk):
        print(self._args)
        return self._func(*[arg._do_calc(chunk) if isinstance(arg, ComputationNode) else arg for arg in self._args])

    def __iter__(self):
        for chunk in self._origin_stream:
            yield self._do_calc(chunk)


class UfuncNode(ComputationNode):
    pass


class ArrayFuncNode(ComputationNode):
    pass 

class LeafNode(ComputationNode):
    def __init__(self, origin_stream, attribute_name):
        self._attribute_name = attribute_name
        self._origin_stream = origin_stream

    def _do_calc(self, chunk):
        return getattr(chunk, self._attribute_name)
    
def read_file(filename, cls):
    origin_stream = open(filename).read_chunks()
    return cls(*[LeafNode(origin_stream, field.name) for field in dataclasses.fields(cls)])

def func(sequence_with_quality):
    (sequence_with_quality.quality + 10).consume()

def consume(node):
    return list(node)
