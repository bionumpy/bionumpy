import pytest
from bionumpy.arithmetics.geometry import Geometry
from bionumpy.datatypes import Bed6
from bionumpy.util.testing import assert_bnpdataclass_equal

@pytest.fixture
def stranded_intervals():
    return Bed6.from_entry_tuples([
        ('chr1', 15, 30, '.', '.', '+'),
        ('chr1', 20, 40, '.', '.', '+'),
        ('chr1', 20, 40, '.', '.', '-')])

@pytest.fixture
def extended_intervals():
    return Bed6.from_entry_tuples([
        ('chr1', 15, 25, '.', '.', '+'),
        ('chr1', 20, 30, '.', '.', '+'),
        ('chr1', 30, 40, '.', '.', '-')])


@pytest.fixture
def geometry():
    return Geometry({"chr1": 100, "chr2": 50})


def test_extend_to_size(geometry, stranded_intervals, extended_intervals):
    extended = geometry.extend_to_size(stranded_intervals, 10)
    assert_bnpdataclass_equal(extended, extended_intervals)
