from .buffers import buffers, data
from bionumpy.io.wig import WigBuffer
import bionumpy as bnp
import pytest


@pytest.fixture
def wig_buffer():
    return bnp.as_encoded_array(buffers['wig'])


@pytest.fixture
def wig_data():
    return data['wig']


def test_start_comment(wig_buffer, wig_data):
    d = WigBuffer.from_raw_buffer(wig_buffer).get_data()
    bnp.util.testing.assert_bnpdataclass_equal(d, wig_data)
