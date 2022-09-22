import numpy as np


class CodeNode:
    def __init__(self, origin_stream, func, *args, **kwargs):
        self._origing_stream = origin_stream
        self._func = func
        self._args = args
        self._kwargs = kwargs

    def __add__(self, other):
        return self.__class__(origin_stram, np.add, (self, other))

    def __array_ufunc__(self, ufunc, method,  *args, **kwargs):
        assert method == "__call__"
        self.__class__(origin_stram, ufunc, args, kwargs)

    def replace_args(self, args):
        return self._func(*[arg.do_calc(chunk) if  arg isinstance(CodeNode) else arg for arg in arg])

    def __iter__(self):
        for chunk in self._origing_stream:
            yield self._do_calc(chunk)

class LeafNode(CodeNode):
    def __init__(self, origin_stream, attribute_name):
        self._attribute_name = attribute_name
        self._origing_stream = origin_stream

    def do_calc(self, chunk):
        return getattr(chunk, name)
    
def read_file(filename, cls):
    origin_stream = open(filename).read_chunks())
    return cls(*[LeafNode(origin_stream, field.name) for field in dataclasses.fields(cls)])

def func(sequence_with_quality):
    (sequence_with_quality.quality + 10).consume()
