import dataclasses
import types
from numbers import Number

import numpy as np

from bionumpy.io.dump_csv import get_column, join_columns


# from bionumpy import EncodedRaggedArray


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


class LazyBNPDataClass:
    pass

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

    @property
    def buffer(self):
        return self._buffer


def create_lazy_class(dataclass):
    field_names = [field.name for field in dataclasses.fields(dataclass)]
    class NewClass(dataclass, LazyBNPDataClass):
        def __init__(self, item_getter, set_values=None):
            self._itemgetter = item_getter
            self._set_values = set_values or {}
            self._computed = False
            self._data = None

        def __len__(self):
            return self._itemgetter.buffer.count_entries()

        def __repr__(self):
            return self[:10].get_data_object().__repr__().replace('with 10 entries', f'with {len(self)} entries')

        def __str__(self):
            return self[:10].get_data_object().__str__().replace('with 10 entries', f'with {len(self)} entries')

        def __getattr__(self, var_name):
            if var_name in self._set_values:
                return self._set_values[var_name]
            if var_name in field_names:
                return self._itemgetter(var_name)
            return super().__getattr__(var_name)

        def __setattr__(self, key, value):
            if key in ['_itemgetter', '_set_values', '_computed', '_data']:
                return super().__setattr__(key, value)
            self._set_values[key] = value

        def __getitem__(self, idx):
            if isinstance(idx, Number):
                idx = [idx]
                return self[idx].get_data_object()[0]
            new_dict = {key: value[idx] for key, value in self._set_values.items()}
            return self.__class__(self._itemgetter[idx], new_dict)

        def __replace__(self, **kwargs):
            new_dict = {key: value for key, value in self._set_values.items()}
            new_dict.update(kwargs)
            return self.__class__(self._itemgetter, new_dict)

        def __iter__(self):
            return iter(self.get_data_object())

        def get_data_object(self):
            if not self._computed:
                self._data = dataclass(*(getattr(self, field.name) for field in dataclasses.fields(dataclass)))
                self._computed = True
            return self._data

        @classmethod
        def empty(cls):
            return dataclass.empty()

        def __array_function__(self, func, types, args, kwargs):
            assert all(issubclass(t, LazyBNPDataClass) for t in types), types
            if func == np.concatenate:
                objects = [a.get_data_object() for a in args[0]]
                args = (objects, ) + args[1:]
                return func(*args, **kwargs)
            return NotImplemented

        def get_buffer(self):
            if not hasattr(self._itemgetter.buffer, 'get_column_range_as_text') or hasattr(self._itemgetter.buffer, 'SKIP_LAZY'):
                return self._itemgetter.buffer.from_data(self.get_data_object())
            columns = []
            for i, field in enumerate(dataclasses.fields(dataclass)):
                if field.name in self._set_values:
                    columns.append(get_column(self._set_values[field.name], field.type))
                else:
                    columns.append(self._itemgetter.buffer.get_column_range_as_text(i, i+1))
            return join_columns(columns, self._itemgetter.buffer.DELIMITER).ravel()


    NewClass.__name__ = dataclass.__name__
    NewClass.__qualname__ = dataclass.__qualname__
    return NewClass
