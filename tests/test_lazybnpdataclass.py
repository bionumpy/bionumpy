import dataclasses

import numpy as np
from numpy.testing import assert_array_equal

import bionumpy as bnp
import pytest

from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.bnpdataclass.lazybnpdataclass import create_lazy_class, ItemGetter
from bionumpy.util.testing import assert_bnpdataclass_equal, assert_encoded_raggedarray_equal, \
    assert_string_array_equal, assert_encoded_raggedarray_like_equal
from .buffers import combos


@bnpdataclass
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

    def validate_if_not(self):
        pass

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

    def validate_if_not(self):
        pass

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


def test_lazy_with_fasta_buffer():
    text, data, bt = combos['fasta']
    b = bt.from_raw_buffer(bnp.as_encoded_array(text))
    lazy = create_lazy_class(bnp.datatypes.SequenceEntry)(ItemGetter(b, bnp.datatypes.SequenceEntry))
    lazy_name = lazy.name
    data_name = [data.name for data in data]
    assert_encoded_raggedarray_like_equal(lazy_name, data_name)
    # assert_encoded_raggedarray_equal(lazy.name, [data.name for data in data])
    lazy.name = bnp.as_encoded_array(['deg'])
    assert_encoded_raggedarray_equal(lazy.name, ['deg'])


def get_pairwise_data(buffer_type):
    text, data, bt = combos[buffer_type]
    b = bt.from_raw_buffer(bnp.as_encoded_array(text))
    cls = bt.dataclass
    lazy = create_lazy_class(bt.dataclass)(ItemGetter(b, bt.dataclass))
    tuples = [tuple(getattr(row, field.name) for field in dataclasses.fields(row)) for row in data]
    return lazy, cls.from_entry_tuples(tuples)


def test_lazy_with_bed_buffer():
    text, data, bt = combos['bed']
    b = bt.from_raw_buffer(bnp.as_encoded_array(text))
    lazy = create_lazy_class(bnp.datatypes.Bed6)(ItemGetter(b, bnp.datatypes.Bed6))
    assert_encoded_raggedarray_like_equal(lazy.chromosome, [data.chromosome for data in data])
    lazy.chromosome = bnp.as_encoded_array(['deg'])
    assert_encoded_raggedarray_equal(lazy.chromosome, ['deg'])


@pytest.fixture
def pair():
    return get_pairwise_data('bed')


@pytest.fixture
def fastq_pair():
    return get_pairwise_data('bed')


@pytest.mark.parametrize('pair', [get_pairwise_data('fasta'), get_pairwise_data('bed')])
@pytest.mark.parametrize('idx', [[0], [1], [0, 1]])
def test_indexing(pair, idx):
    lazy, dataclass = pair
    assert_bnpdataclass_equal(lazy.get_data_object(), dataclass)
    data_object = lazy[idx].get_data_object()
    assert_bnpdataclass_equal(data_object, dataclass[idx])


@pytest.mark.parametrize('idx', [[0], [1], [0, 1]])
def test_mutated_indexing(idx):
    lazy, dataclass = get_pairwise_data('bed')
    new_starts = np.full(len(dataclass), 100000)
    lazy.start = new_starts
    dataclass.start = new_starts
    assert_bnpdataclass_equal(lazy.get_data_object(), dataclass)
    data_object = lazy[idx].get_data_object()
    assert_bnpdataclass_equal(data_object, dataclass[idx])

@pytest.mark.parametrize('pair', [get_pairwise_data('fasta'), get_pairwise_data('bed')])
def test_str(pair):
    lazy, dataclass = pair
    assert str(lazy) == str(dataclass)
    assert repr(lazy) == str(dataclass)
