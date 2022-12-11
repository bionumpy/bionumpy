from bionumpy.util.ascii_hash import get_ascii_hash, column_index_array, AsciiHashTable
import pytest
import bionumpy as bnp
from numpy.testing import assert_array_equal
from npstructures.testing import assert_raggedarray_equal


@pytest.fixture
def encoded_array():
    return bnp.as_encoded_array("chr1")


@pytest.fixture
def encoded_ragged_array():
    return bnp.as_encoded_array(["chr1",
                                 "chr2",
                                 "chr13"])


def test_ascii_hash_on_encoded_array(encoded_array):
    """
    In [5]: (99+104*128+114*128**2+49*128**3) % 103
    Out[5]: 21
    """
    h = get_ascii_hash(encoded_array, 103)
    assert h == 21


def test_column_index_array(encoded_ragged_array):
    assert_raggedarray_equal(
        column_index_array(encoded_ragged_array.shape),
        [[0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3, 4]])


def test_ascii_hash_on_encoded_ragged_array(encoded_ragged_array):
    """
    In [5]: (99+104*128+114*128**2+49*128**3) % 103
    Out[5]: 21
    """
    h = get_ascii_hash(encoded_ragged_array, 103)
    assert_array_equal(h, [21, 93, 48])


def test_ascii_string_hash_table(encoded_ragged_array):
    hash_table = AsciiHashTable.from_sequences(encoded_ragged_array)
    assert_array_equal(hash_table[encoded_ragged_array],
                       [0, 1, 2])
    
