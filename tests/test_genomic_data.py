import numpy as np
from numpy.testing import assert_equal
from bionumpy.datatypes import Interval, StrandedInterval
from bionumpy.arithmetics.intervals import GenomicRunLengthArray
from bionumpy.genomic_data.genomic_intervals import GenomicIntervals, GenomicLocation
from bionumpy.genomic_data.genomic_track import GenomicArray, GenomicArrayGlobal, GenomicArrayNode
from bionumpy.util.testing import assert_raggedarray_equal
from bionumpy.computation_graph import compute
import pytest


@pytest.fixture
def chrom_sizes():
    return {"chr1": 100, "chr2": 100}


@pytest.fixture
def intervals():
    return Interval.from_entry_tuples([
        ('chr1', 10, 20),
        ('chr2', 20, 30),
        ('chr2', 50, 70)])


@pytest.fixture
def stranded_intervals():
    return StrandedInterval.from_entry_tuples([
        ('chr1', 10, 20, '-'),
        ('chr2', 20, 30, '+'),
        ('chr2', 50, 70, '+')])


@pytest.fixture
def windows(chrom_sizes):
    w =  GenomicIntervals.from_fields(
        chrom_sizes,
        ['chr1', 'chr2', 'chr2'],
        [14, 15, 45],
        [24, 25, 55],
        ['-', '+', '+'])
    assert w.is_stranded()
    return w


@pytest.fixture
def genomic_locations(chrom_sizes):
    return GenomicLocation.from_fields(
        chrom_sizes,
        ['chr1', 'chr2', 'chr2'],
        [19, 20, 50],
        ['-', '+', '+'])


@pytest.fixture
def track_arrays(chrom_sizes):
    return {name: np.repeat(np.arange(20), 5)[:size]
            for name, size in chrom_sizes.items()}


@pytest.fixture
def track(track_arrays):
    return GenomicArrayGlobal.from_dict(
        {name: GenomicRunLengthArray.from_array(a)
         for name, a in track_arrays.items()})


@pytest.fixture
def track_stream(track_arrays):
    return GenomicArrayNode.from_dict(
        {name: GenomicRunLengthArray.from_array(a)
         for name, a in track_arrays.items()})


@pytest.fixture
def genomic_intervals(stranded_intervals, chrom_sizes):
    return GenomicIntervals.from_intervals(stranded_intervals, chrom_sizes, is_stranded=True)


def test_get_location(genomic_intervals, genomic_locations):
    result = genomic_intervals.get_location('start')
    assert_equal(result.position, genomic_locations.position)


def test_get_windows(genomic_locations, windows):
    result = genomic_locations.get_windows(flank=5)
    assert_equal(result.start, windows.start)
    assert_equal(result.stop, windows.stop)


def test_extract_windows(windows, track_stream, track_arrays):
    tmp = track_stream[windows].mean(axis=0)
    result, = compute(tmp)
    result = result.to_array()
    true = np.mean([track_arrays[w.chromosome.to_string()][w.start:w.stop] if w.strand == '+' else track_arrays[w.chromosome.to_string()][w.stop-1:w.start-1:-1] for w in windows], axis=0)
    assert_raggedarray_equal(result, true)
                             
