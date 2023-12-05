from numpy.testing import assert_array_equal

import bionumpy as bnp
from bionumpy.genomic_data.coordinate_mapping import find_indices


def test_find_indices():
    locations = [2, 3, 5, 7, 11]
    intervals = bnp.Interval.from_entry_tuples(
        [('chr1', 1, 8), ('chr1', 3, 6), ('chr1', 5, 10)])
    location_indices, interval_indices = find_indices(locations, intervals)
    true_interval_indices = [0, 0, 0, 0,
                             1, 1,
                             2, 2]
    true_location_indices = [0, 1, 2, 3,
                             1, 2,
                             2, 3]
    assert_array_equal(interval_indices, true_interval_indices)
    assert_array_equal(location_indices, true_location_indices)
