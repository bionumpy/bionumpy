from bionumpy.streams import MultiStream
from bionumpy.arithmetics.similarity_measures import forbes, jaccard, get_contingency_table
from bionumpy.datatypes import Interval
from numpy.testing import assert_array_equal
import pytest


@pytest.fixture
def interval_a():
    return Interval.from_entry_tuples([("chr1", 10, 20), ("chr2", 20, 30)])


@pytest.fixture
def interval_b():
    return Interval.from_entry_tuples([("chr1", 15, 22), ("chr2", 15, 25)])


@pytest.fixture
def chromosome_sizes():
    return {"chr1": 100, "chr2": 50}


def test_contigency_table(interval_a, interval_b, chromosome_sizes):
    ms = MultiStream(chromosome_sizes, a=interval_a, b=interval_b)
    contigency_table = get_contingency_table(ms.a, ms.b, ms.lengths)
    assert_array_equal(contigency_table, [[10, 10], [7, 123]])


def test_forbes(interval_a, interval_b, chromosome_sizes):
    f = forbes(chromosome_sizes, interval_a, interval_b)
    true = (150*10)/(20*17)
    assert f == true


def test_jaccard(interval_a, interval_b, chromosome_sizes):
    f = jaccard(chromosome_sizes, interval_a, interval_b)
    true = 10/(12+15)
    assert f == true
