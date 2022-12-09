from bionumpy.util.ascii_hash import get_ascii_hash
import pytest
import bionumpy as bnp


@pytest.fixture
def encoded_array():
    return bnp.as_encoded_array("chr1")


def test_ascii_hash_on_encoded_array(encoded_array):
    """
    In [5]: (99+104*128+114*128**2+49*128**3) % 103
    Out[5]: 21
    """
    h = get_ascii_hash(encoded_array, 103)
    assert h == 21
