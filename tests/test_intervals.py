import numpy as np
import pytest
from bionumpy.arithmetics import count_overlap, intersect, get_pileup, sort_intervals, get_boolean_mask, unique_intersect
from bionumpy.util.testing import assert_bnpdataclass_equal
from bionumpy.datatypes import Interval, BedGraph


@pytest.fixture
def interval_a():
    return Interval(["chr1"]*3, [10, 20, 30], [15, 29, 35])


@pytest.fixture
def small_interval():
    return Interval.from_entry_tuples([
        ("chr1", 2, 5),
        ("chr1", 7, 9)])


@pytest.fixture
def small_complicated_interval():
    return Interval.from_entry_tuples([
        ("chr1", 2, 5),
        ("chr1", 5, 7),
        ("chr1", 10, 12),
        ("chr1", 11, 13)])


def complicated_intervals():
    return [Interval.from_entry_tuples([
        ("chr1", 2, 5),
        ("chr1", 5, 7),
        ("chr1", 10, 12),
        ("chr1", 11, 13)]),
     Interval.from_entry_tuples([
         ("chr1", 0, 5),
         ("chr1", 11, 13)]),
     Interval.from_entry_tuples([
         ("chr1", 0, 5),
         ("chr1", 11, 20)])]


@pytest.fixture
def interval_b():
    return Interval(["chr1"]*3, [10, 22, 29], [15, 28, 36])


@pytest.fixture
def interval_c():
    return Interval(["chr1"]*3, [10, 15, 29], [15, 28, 36])


@pytest.fixture
def interval_d():
    return Interval(["chr3", "chr2", "chr2", "chr1"],
                    [10, 15, 14, 12],
                    [20, 22, 23, 24])


def test_count_overlap(interval_a, interval_b):
    assert count_overlap(interval_a, interval_b) == 5+6+5


def test_intersect(interval_a, interval_b):
    true = Interval(["chr1"]*3, [10, 22, 30], [15, 28, 35])
    result = intersect(interval_a, interval_b)
    assert_bnpdataclass_equal(true, result)

@pytest.fixture
def small_complicated_interval():
    return Interval.from_entry_tuples([
        ("chr1", 2, 5),
        ("chr1", 5, 7),
        ("chr1", 10, 12),
        ("chr1", 11, 13)])


def test_unique_intersect(small_complicated_interval):
    ref = Interval.from_entry_tuples([('chr1', 7, 11)])
    true = small_complicated_interval[[2]]
    result = unique_intersect(small_complicated_interval, ref, 20)
    print(result)
    assert_bnpdataclass_equal(true, result)

@pytest.mark.skip("obsolete")
def test_pileup(interval_c):
    p = pileup(interval_c)
    print("#", p.chromosome)
    assert_bnpdataclass_equal(p == BedGraph(["chr1"]*3, [10, 28, 29], [28, 29, 36], [1, 0, 1]))


def test_sort_intervals(interval_d):
    sorted_intervals = sort_intervals(interval_d)
    assert_bnpdataclass_equal(sorted_intervals,
                              interval_d[[3, 2, 1, 0]])


@pytest.mark.parametrize("interval", complicated_intervals())
def test_get_boolean_mask(interval):
    size = 20
    starts = np.asanyarray(interval.start)
    ends = np.asanyarray(interval.stop)
    row_len = size
    true = np.zeros((row_len), dtype=bool)
    for (start, end) in zip(starts, ends):
        print(start, end)
        true[start:end] |= True
    tmp = get_boolean_mask(interval, size)
    print(tmp.starts, tmp.ends, tmp.values)
    print(tmp._events, tmp._values)
    result = tmp.to_array()
    np.testing.assert_array_equal(result, true)

# @pytest.mark.parametrize("interval", complicated_intervals())
# def test_merge_intervals(interval):
#     size = 20
#     starts = np.asanyarray(interval.start)
#     ends = np.asanyarray(interval.stop)
#     row_len = size
#     true = np.zeros((row_len), dtype=bool)
#     for (start, end) in zip(starts, ends):
#         print(start, end)
#         true[start:end] |= True
#     tmp = get_boolean_mask(interval, size)
#     print(tmp.starts, tmp.ends, tmp.values)
#     result = tmp.to_array()
#     np.testing.assert_array_equal(result, true)

