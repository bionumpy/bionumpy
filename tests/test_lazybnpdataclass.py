import dataclasses

import numpy as np
from numpy.testing import assert_array_equal

import bionumpy as bnp
import pytest

from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.bnpdataclass.lazybnpdataclass import create_lazy_class, ItemGetter
from bionumpy.util.testing import assert_bnpdataclass_equal, assert_encoded_raggedarray_equal


@dataclasses.dataclass
class DummyClass:
    a: int
    b: str


@bnpdataclass
class DummyBNP:
    a: int
    b: str

@pytest.mark.skip
def test_buffer_backed_descriptor():
    buffer = [0, 1, 2]

    class NewClass:
        attr_a: BufferBackedDescriptor(buffer, 0, str)
        attr_b: BufferBackedDescriptor(buffer, 1, int)
        attr_c: BufferBackedDescriptor(buffer, 2, float)


def itemgetter(dataclass):
    str_getter = lambda i: f'str({i})'
    int_getter = lambda i: i * 2
    return {str: str_getter, int: int_getter}


def silly_property_getter():
    i = 0

    def my_prop():
        nonlocal i
        i += 1
        return i

    return my_prop()


def simple_backend(var_name):
    if var_name == 'a':
        return 10
    elif var_name == 'b':
        return 'hei'
    return None


class BufferMock:
    def get_field_by_number(self, i: int, t: type):
        if i == 0:
            return 10
        elif i == 1:
            return 'hei'
        return None


class BetterBufferMock:
    def __init__(self, *data):
        self.data = data

    def get_field_by_number(self, i, t):
        return self.data[i]


@pytest.mark.parametrize('backend', [simple_backend, ItemGetter(BufferMock(), DummyClass)])
def test_create_lazy_class(backend):
    lazy_obj = create_lazy_class(DummyClass)(backend)
    assert lazy_obj.a == 10
    assert lazy_obj.b == 'hei'
    lazy_obj.b = 20
    assert lazy_obj.b == 20


@pytest.fixture
def lazy_bnp():
    b = bnp.as_encoded_array(['hei'])
    a = np.array([10])
    mock = BetterBufferMock(a, b)
    return create_lazy_class(DummyBNP)(ItemGetter(mock, DummyBNP))


def test_bnp_lazy(lazy_bnp):
    assert_array_equal(lazy_bnp.a, [10])
    assert_encoded_raggedarray_equal(lazy_bnp.b, ['hei'])
    lazy_bnp.a = [20]
    assert_array_equal(lazy_bnp.a, [20])
