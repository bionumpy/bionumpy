import pytest
from bionumpy.intervals import count_overlap, intersect
from bionumpy.datatypes import Interval


@pytest.fixture
def interval_a():
    return Interval(["chr1"]*3, [10, 20, 30], [15, 29, 35])


@pytest.fixture
def interval_b():
    return Interval(["chr1"]*3, [10, 22, 29], [15, 28, 36])


def test_count_overlap(interval_a, interval_b):
    assert count_overlap(interval_a, interval_b) == 5+6+5


def test_intersect(interval_a, interval_b):
    true = Interval(["chr1"]*3, [10, 22, 30], [15, 28, 35])
    result = intersect(interval_a, interval_b)
    print(result)
    print(true)
    assert true == result
