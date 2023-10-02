from bionumpy.datatypes import Interval
from bionumpy.arithmetics import get_pileup
from numpy.testing import assert_array_equal
import numpy as np
import pytest

size = 12


interval_list = [
    Interval(["chr1"]*4, [2, 3, 5, 7], [4, 6, 8, 10]),
    Interval(["chr1"]*5, [2, 3, 3, 5, 7], [4, 6, 8, 8, 10]),
    Interval(["chr1"]*5, [2, 3, 3, 5, 7], [4, 6, 7, 8, 10])]


def raw_pileup(intervals, size):
    pileup = np.zeros(size, dtype=int)
    for interval in intervals:
        pileup[interval.start:interval.stop] += 1
    return pileup


@pytest.mark.parametrize("intervals", interval_list)
def test_pileup(intervals):
    p = get_pileup(intervals, size)
    raw_p = raw_pileup(intervals, size)
    assert_array_equal(p.to_array(), raw_p)
