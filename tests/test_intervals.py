import pytest
from bionumpy.intervals import count_overlap, intersect, pileup
from bionumpy.datatypes import Interval
from bionumpy.bedgraph import BedGraph

@pytest.fixture
def interval_a():
    return Interval(["chr1"]*3, [10, 20, 30], [15, 29, 35])


@pytest.fixture
def interval_b():
    return Interval(["chr1"]*3, [10, 22, 29], [15, 28, 36])

@pytest.fixture
def interval_c():
    return Interval(["chr1"]*3, [10, 15, 29], [15, 28, 36])


def test_count_overlap(interval_a, interval_b):
    assert count_overlap(interval_a, interval_b) == 5+6+5


def test_intersect(interval_a, interval_b):
    true = Interval(["chr1"]*3, [10, 22, 30], [15, 28, 35])
    result = intersect(interval_a, interval_b)
    assert true == result

@pytest.mark.skip("obsolete")
def test_pileup(interval_c):
    p = pileup(interval_c)
    print("#", p.chromosome)
    assert p == BedGraph(["chr1"]*3, [10, 28, 29], [28, 29, 36], [1, 0, 1])
