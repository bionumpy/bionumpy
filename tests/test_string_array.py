import numpy as np
import pytest
from numpy.testing import assert_array_equal
import bionumpy as bnp
from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.io.delimited_buffers import DelimitedBuffer
from bionumpy.typing import SequenceID
from bionumpy.string_array import StringArray, string_array


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
def file_name():
    name = 'string_array_test.txt'
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
