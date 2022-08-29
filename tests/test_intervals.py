import pytest
from bionumpy.intervals import count_overlap
from bionumpy.datatypes import Interval


@pytest.fixture
def interval_a():
    return Interval(None, [10, 20, 30], [15, 29, 35])


@pytest.fixture
def interval_b():
    return Interval(None, [10, 22, 29], [15, 28, 36])


def test_count_overlap(interval_a, interval_b):
    assert count_overlap(interval_a, interval_b) == 5+6+5
