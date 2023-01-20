import pytest
import numpy as np
from bionumpy import str_equal
from bionumpy.streams import groupby
from bionumpy.datatypes import Interval, BedGraph
from numpy.testing import assert_equal
from bionumpy.arithmetics.geometry import Geometry, GenomicTrack, StreamedGeometry
from bionumpy.datatypes import Bed6, ChromosomeSize
from bionumpy.util.testing import assert_bnpdataclass_equal


@pytest.fixture
def stranded_intervals():
    return Bed6.from_entry_tuples([
        ('chr1', 15, 30, '.', '.', '+'),
        ('chr1', 20, 40, '.', '.', '+'),
        ('chr1', 20, 40, '.', '.', '-'),
        ('chr2', 10, 20, '.', '.', '+')])


@pytest.fixture
def extended_intervals():
    return Bed6.from_entry_tuples([
        ('chr1', 15, 25, '.', '.', '+'),
        ('chr1', 20, 30, '.', '.', '+'),
        ('chr1', 30, 40, '.', '.', '-'),
        ('chr2', 10, 20, '.', '.', '+')])


@pytest.fixture
def invalid_intervals():
    return Bed6.from_entry_tuples([
        ('chr1', -10, 10, '.', '.', '+'),
        ('chr1', 90, 110, '.', '.', '+'),
        ('chr1', 30, 40, '.', '.', '-'),
        ('chr2', -10, 60, '.', '.', '+')])


@pytest.fixture
def valid_intervals():
    return Bed6.from_entry_tuples([
        ('chr1', 0, 10, '.', '.', '+'),
        ('chr1', 90, 100, '.', '.', '+'),
        ('chr1', 30, 40, '.', '.', '-'),
        ('chr2', 0, 50, '.', '.', '+')])


@pytest.fixture
def disjoint_intervals():
    return Bed6.from_entry_tuples([
        ('chr1', 0, 10, '.', '.', '+'),
        ('chr1', 90, 100, '.', '.', '+'),
        ('chr2', 0, 50, '.', '.', '+')])


@pytest.fixture
def pileup():
    return BedGraph.from_entry_tuples([
        ('chr1', 0, 10, 1),
        ('chr1', 10, 90, 2),
        ('chr1', 90, 100, 3),
        ('chr2', 0, 50, 4)])


@pytest.fixture
def chrom_sizes():
    return {"chr1": 100, "chr2": 50}

@pytest.fixture
def geometry(chrom_sizes):
    return Geometry(chrom_sizes)

@pytest.fixture
def streamed_geometry(chrom_sizes):
    return StreamedGeometry(chrom_sizes)


def test_from_chrom_sizes(chrom_sizes):
    c = ChromosomeSize.from_entry_tuples(
        list(chrom_sizes.items()))
    g = Geometry.from_chrom_sizes(c)
    assert g.chrom_size('chr2') == 50

def test_clip(geometry, invalid_intervals, valid_intervals):
    clipped = geometry.clip(invalid_intervals)
    assert_bnpdataclass_equal(clipped, valid_intervals)


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


def test_to_intervals(geometry, disjoint_intervals):
    mask = geometry.get_mask(disjoint_intervals)
    intervals = mask.to_intervals()
    assert_bnpdataclass_equal(intervals, disjoint_intervals.astype(Interval))


def test_to_bedgraph(geometry, pileup):
    track = geometry.get_track(pileup)
    new_pileup = track.to_bedgraph()
    assert_bnpdataclass_equal(pileup, new_pileup)


def test_stream_sum(pileup, geometry, streamed_geometry):
    track = streamed_geometry.get_track(pileup)
    assert track.sum() == geometry.get_track(pileup).sum()


def test_compute(pileup, geometry, streamed_geometry):
    track = streamed_geometry.get_track(pileup).compute()
    assert_equal(track, geometry.get_track(pileup))
