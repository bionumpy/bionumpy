import pytest
from npstructures.testing import assert_raggedarray_equal
import numpy as np
from bionumpy import as_sequence_array
from bionumpy.strops import (int_to_str, ints_to_strings, join, split, str_to_int, str_equal)


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


def test_str_to_int(ints, strings):
    for i, s in zip(ints, str_to_int(strings)):
        assert i == s


def test_join(strings):
    seqs = as_sequence_array(strings)
    joined = join(seqs)
    true = as_sequence_array("\t".join(strings))
    assert np.all(joined == as_sequence_array("\t".join(strings)))


def test_split(strings):
    joined = as_sequence_array(",".join(strings))
    assert_raggedarray_equal(split(joined), as_sequence_array(strings))


def test_str_equal(strings):
    s = as_sequence_array(strings*2)
    mask = str_equal(s, "12")
    np.testing.assert_array_equal(mask, [False, True, False]*2)
