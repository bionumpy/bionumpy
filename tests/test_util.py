import pytest
import numpy as np
#from bionumpy.util import filter_on_intervals
from bionumpy.arithmetics import sort_intervals, merge_intervals
from npstructures import npdataclass


@npdataclass
class Entry:
    position: np.ndarray


@npdataclass
class Interval:
    start: np.ndarray
    stop: np.ndarray


@pytest.fixture
def entry():
    return Entry(np.arange(13))


@pytest.fixture
def intervals():
    return Interval(np.array([2, 7, 11]), np.array([5, 11, 13]))


@pytest.fixture
def sorted_intervals():
    return Interval(np.array([1, 2, 2, 10, 12, 20]),
                    np.array([5, 3, 4, 15, 17, 25]))

@pytest.fixture
def unsorted_intervals():
    return Interval(np.array([1, 2, 10, 2, 12, 20][::-1]),
                    np.array([5, 4, 15, 3, 17, 25][::-1]))


@pytest.fixture
def merged_intervals():
    return Interval(np.array([1, 10, 20]),
                    np.array([5, 17, 25]))

@pytest.mark.skip("outdated")
def test_sort_intervals(unsorted_intervals, sorted_intervals):
    s = sort_intervals(unsorted_intervals)
    assert s == sorted_intervals


def test_merged_intervals(sorted_intervals, merged_intervals):
    m = merge_intervals(sorted_intervals)
    assert m == merged_intervals


@pytest.mark.skip
def test_filter_on_intervals(entry, intervals):
    truth = Entry([2, 3, 4, 7, 8, 9, 10, 11, 12])
    assert np.all(filter_on_intervals(entry, intervals) == truth)
