import pytest
from npstructures.testing import assert_raggedarray_equal
import numpy as np
from bionumpy import as_sequence_array
from bionumpy.strops import int_to_str, ints_to_strings, join


@pytest.fixture()
def ints():
    return [1, 12, 123]


@pytest.fixture()
def strings():
    return ["1", "12", "123"]


def test_int_to_string(ints, strings):
    for i, s in zip(ints, strings):
        assert str(int_to_str(i)) == s


def test_multiple_int_to_string(ints, strings):
    for i, s in zip(ints_to_strings(ints), strings):
        assert str(i) == s


def test_join(strings):
    seqs = as_sequence_array(strings)
    joined = join(seqs)
    true = as_sequence_array("\t".join(strings))
    assert np.all(joined == as_sequence_array("\t".join(strings)))
