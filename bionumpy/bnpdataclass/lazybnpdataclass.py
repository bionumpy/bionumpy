import dataclasses
import types

import numpy as np

from bionumpy import EncodedRaggedArray


def translate_types(input_type):
    if input_type == str:
        return EncodedRaggedArray
    elif input_type == int:
        return np.ndarray


def buffer_backed_bnp(old_cls):
    cls = types.new_class(old_cls.__name__, old_cls.__bases__, {})
    for i, (var_name, var_type) in enumerate(old_cls.__annotations__.items()):
        setattr(cls, var_name, BufferBackedDescriptor(i, var_type))
        setattr(cls, '__init__', lambda self, buffer: setattr(self, '_buffer', buffer))
    return cls


class BufferBackedDescriptor:
    '''
    This class is made to access and parse parts of a text buffer lazily.v
    '''

    def __init__(self, buffer, index, dtype):
        self._buffer = buffer
        self._index = index
        self._dtype = dtype

    def __get__(self, obj, objtype):
        return self._dtype(obj._buffer.get_field_by_number(self._index, self._dtype))


class BaseClass:
    def __init__(self, buffer):
        self._buffer = buffer

    def __getattr__(self, var_name):
        if var_name in self._buffer:
            return self._buffer[var_name]
        return super().__getattr__(var_name)


class ItemGetter:
    def __init__(self, buffer: 'FileBuffer', dataclass: dataclasses.dataclass):
        self._buffer = buffer
        self._dataclass = dataclass
        self._field_dict = {field.name: (i, field.type) for i, field in enumerate(dataclasses.fields(dataclass))}

    def __call__(self, name):
        return self._buffer.get_field_by_number(self._field_dict[name][0], self._field_dict[name][1])

    def __getitem__(self, idx):
        return self.__class__(self._buffer[idx], self._dataclass)


def create_lazy_class(dataclass):
    class NewClass(dataclass):
        def __init__(self, item_getter, set_values=None):
            self._itemgetter = item_getter
            self._set_values = set_values or {}

        def __repr__(self):
            return self[:10].get_data_object().__repr__()

        def __str__(self):
            return self[:10].get_data_object().__str__()

        def __getattr__(self, var_name):
            if var_name in self._set_values:
                return self._set_values[var_name]
            return self._itemgetter(var_name)

        def __setattr__(self, key, value):
            if key in ['_itemgetter', '_set_values']:
                return super().__setattr__(key, value)
            self._set_values[key] = value

        def __getitem__(self, idx):
            new_dict = {key: value[idx] for key, value in self._set_values.items()}
            return self.__class__(self._itemgetter[idx], new_dict)

        def get_data_object(self):
            return dataclass(*(getattr(self, field.name) for field in dataclasses.fields(dataclass)))

    NewClass.__name__ = dataclass.__name__
    NewClass.__qualname__ = dataclass.__qualname__
    return NewClass
