import dataclasses
import types
from functools import lru_cache
from numbers import Number

import numpy as np

from bionumpy.io.dump_csv import get_column, join_columns
from bionumpy.io.exceptions import FormatException


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
    This class is made to access and parse parts of a text buffer lazily.
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
    def __init__(self, buffer: 'FileBuffer', dataclass: dataclasses.dataclass, start_line=0):
        self._buffer = buffer
        self._dataclass = dataclass
        self._field_dict = {field.name: (i, field.type)
                            for i, field in
                            enumerate(dataclasses.fields(dataclass))}
        self._buffer.validate_if_not()
        self._start_line = start_line

    def concatenate(self, itemgetters):
        return self.__class__(self._buffer.concatenate([i._buffer for i in itemgetters]),
                              itemgetters[0]._dataclass, itemgetters[0]._start_line)

    @lru_cache()
    def n_entries(self):
        return self._buffer.count_entries()

    def __call__(self, name):
        try:
            return self._buffer.get_field_by_number(self._field_dict[name][0], self._field_dict[name][1])
        except FormatException as e:
            e.line_number += self._start_line
            raise e

    def __getitem__(self, idx):
        return self.__class__(self._buffer[idx], self._dataclass)

    @property
    def buffer(self):
        return self._buffer


def create_lazy_class(dataclass, header=None):
    field_names = [field.name for field in dataclasses.fields(dataclass)]

    class NewClass(dataclass, LazyBNPDataClass):
        def __init__(self, item_getter, set_values=None, computed_values=None):
            self._itemgetter = item_getter
            self._set_values = set_values or {}
            self._computed_values = computed_values or {}
            self._computed = False
            self._data = None
            self._header = header
            # if header is not None:
            # self.set_context('header', header)

        @classmethod
        def from_data_frame(cls, df):
            return dataclass.from_data_frame(df)

        @classmethod
        def from_dict(cls, d):
            return dataclass.from_dict(d)

        def toiter(self):
            return self.get_data_object().toiter()

        def tolist(self):
            return self.get_data_object().tolist()

        def todict(self):
            return self.get_data_object().todict()

        def topandas(self):
            return self.get_data_object().topandas()

        def __len__(self):
            return self._itemgetter.n_entries()

        def __repr__(self):
            return self[:10].get_data_object().__repr__().replace('with 10 entries', f'with {len(self)} entries')

        def __str__(self):
            return self[:10].get_data_object().__str__().replace('with 10 entries', f'with {len(self)} entries')

        def __getattr__(self, var_name):
            if var_name in self._set_values:
                return self._set_values[var_name]
            if var_name in field_names:
                if var_name not in self._computed_values:
                    value = self._itemgetter(var_name)
                    self._computed_values[var_name] = value
                return self._computed_values[var_name]
            return getattr(super(), var_name)
            # raise ValueError(f'No such field {var_name} in {self.__class__.__name__}')

        def __setattr__(self, key, value):
            if key in ['_itemgetter', '_set_values', '_computed', '_data', '_computed_values', '_header']:
                return super().__setattr__(key, value)
            self._set_values[key] = value
            if key in self._computed_values:
                del self._computed_values[key]

        def __getitem__(self, idx):
            if isinstance(idx, Number):
                idx = [idx]
                return self[idx].get_data_object()[0]
            new_dict = {key: value[idx] for key, value in self._set_values.items()}
            new_computed = {key: value[idx] for key, value in self._computed_values.items()}
            return self.__class__(self._itemgetter[idx], new_dict, new_computed)

        def __replace__(self, **kwargs):
            new_dict = {key: value for key, value in self._set_values.items()}
            new_dict.update(kwargs)
            return self.__class__(self._itemgetter, new_dict)

        def __iter__(self):
            return iter(self.get_data_object())

        def get_data_object(self):
            if not self._computed:
                fields = [getattr(self, field.name) for field in dataclasses.fields(dataclass)]
                self._data = dataclass(*fields)
                self._computed = True
            return self._data

        @classmethod
        def empty(cls):
            return dataclass.empty()

        def __array_function__(self, func, types, args, kwargs):
            assert all(issubclass(t, LazyBNPDataClass) for t in types), types
            if func == np.concatenate:
                values = args[0]
                if hasattr(values[0]._itemgetter.buffer, 'concatenate'):
                    if all(not a._set_values for a in values):
                        return self.__class__(self._itemgetter.concatenate([a._itemgetter for a in values]))

                objects = [a.get_data_object() for a in args[0]]
                args = (objects,) + args[1:]
                return func(*args, **kwargs)
            return NotImplemented

        def get_buffer(self, buffer_class=None):
            if buffer_class is None:
                buffer_class = self._itemgetter.buffer.__class__
            if not hasattr(self._itemgetter.buffer, 'get_field_range_as_text') or hasattr(self._itemgetter.buffer,
                                                                                          'SKIP_LAZY') or hasattr(
                    buffer_class, 'SKIP_LAZY'):
                return self._itemgetter.buffer.from_data(self.get_data_object())
            columns = []
            if not self._set_values and issubclass(self._itemgetter.buffer.__class__, buffer_class):
                return self._itemgetter.buffer.data.ravel()
            for i, field in enumerate(dataclasses.fields(dataclass)):
                if field.name in self._set_values:
                    columns.append(get_column(self._set_values[field.name], field.type))
                else:
                    columns.append(self._itemgetter.buffer.get_field_range_as_text(i, i + 1))
            return buffer_class.join_fields(columns)

        def get_context(self, name):
            if name == 'header':
                return self._header

        def has_context(self, name):
            return name == 'header'

    NewClass.__name__ = dataclass.__name__
    NewClass.__qualname__ = dataclass.__qualname__
    return NewClass
