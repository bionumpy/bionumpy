import numpy as np
import pytest
from numpy.testing import assert_array_equal
import bionumpy as bnp
from bionumpy import as_encoded_array
from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.io.delimited_buffers import DelimitedBuffer
from bionumpy.typing import SequenceID
from bionumpy.string_array import StringArray, string_array
from bionumpy.util.testing import assert_bnpdataclass_equal


@bnpdataclass
class SimpleEntry:
    name: SequenceID
    age: int

class SimpleBuffer(DelimitedBuffer):
    dataclass = SimpleEntry


def test_from_strings():
    input_strings = ['hei', 'pa', 'deg']
    string_array = StringArray(input_strings)
    assert_array_equal(string_array.raw(), np.array([b'hei', b'pa', b'deg'], dtype='S'))
    assert np.all(string_array == input_strings)


def test_bnpclass():
    se = SimpleEntry(name=['hei', 'pa', 'deg'], age=[1, 2, 3])
    raw = se.name.raw()
    assert_array_equal(raw, np.array([b'hei', b'pa', b'deg'], dtype='S'))


@pytest.fixture
def file_name(tmp_path):
    name = tmp_path / 'string_array_test.txt'
    open(name, 'w').write(
        '''\
hei\t1
pa\t2
deg\t3
''')
    return name


def test_read(file_name):
    data = bnp.open(file_name, buffer_type=SimpleBuffer).read()
    assert np.all(data.name == ['hei', 'pa', 'deg'])
    data = data.get_data_object()
    with bnp.open(file_name, 'w', buffer_type=SimpleBuffer) as f:
        f.write(data)
    new_data = bnp.open(file_name, buffer_type=SimpleBuffer).read()
    assert_bnpdataclass_equal(data, new_data)

def test_isin():
    array = string_array(['hei', 'pa', 'deg'])
    assert_array_equal(np.isin(array, ['hei', 'pa']), np.array([True, True, False]))

def test_empty():
    a = as_encoded_array([])
    s = string_array(a)
    assert len(s) == 0