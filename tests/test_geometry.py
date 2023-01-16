import pytest
import numpy as np
from bionumpy import str_equal
from numpy.testing import assert_equal
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


def test_get_pileup(geometry, stranded_intervals):
    genomic_track = geometry.get_pileup(stranded_intervals)
    for chromosome, track in genomic_track.to_dict().items():
        true = np.zeros(geometry.chrom_size(chromosome), dtype=int)
        for interval in stranded_intervals:
            if str_equal(interval.chromosome, chromosome):
                true[interval.start:interval.stop] += 1
        assert_equal(true, track)
