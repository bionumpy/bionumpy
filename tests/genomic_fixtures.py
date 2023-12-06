import numpy as np
from numpy.testing import assert_equal
from bionumpy.datatypes import Interval, StrandedInterval, BedGraph
from bionumpy.arithmetics.intervals import GenomicRunLengthArray
from bionumpy.genomic_data.genomic_intervals import GenomicIntervals, GenomicLocation, LocationEntry
from bionumpy.genomic_data.genomic_track import GenomicArray, GenomicArrayGlobal, GenomicArrayNode
from bionumpy.genomic_data.genome_context import GenomeContext, GenomeError
from bionumpy.genomic_data.genome import Genome
from bionumpy.util.testing import assert_raggedarray_equal, assert_bnpdataclass_equal
from bionumpy.computation_graph import compute
import pytest


@pytest.fixture
def chrom_sizes():
    return {"chr1": 100, "chr2": 100}


@pytest.fixture
def genome_context(chrom_sizes):
    return GenomeContext.from_dict(chrom_sizes)


@pytest.fixture
def chrom_sizes_big():
    return {'chr1': 10,
            'chr2': 20,
            'chr10_alt': 100,
            'chr3': 30}


@pytest.fixture
def genome_context_big(chrom_sizes_big):
    return GenomeContext.from_dict(chrom_sizes_big)


@pytest.fixture
def intervals():
    return Interval.from_entry_tuples([
        ('chr1', 10, 20),
        ('chr2', 20, 30),
        ('chr2', 50, 70)])


@pytest.fixture
def intervals_big():
    return Interval.from_entry_tuples([
        ('chr1', 5, 7),
        ('chr10_alt', 1, 2),
        ('chr3', 2, 3),
        ('chr3', 5, 7)])


@pytest.fixture
def stranded_intervals():
    return StrandedInterval.from_entry_tuples([
        ('chr1', 10, 20, '-'),
        ('chr2', 20, 30, '+'),
        ('chr2', 50, 70, '+')])


@pytest.fixture
def windows(genome_context):
    w = GenomicIntervals.from_fields(
        genome_context,
        ['chr1', 'chr2', 'chr2'],
        [14, 15, 45],
        [25, 26, 56],
        ['-', '+', '+'])
    assert w.is_stranded()
    return w


@pytest.fixture
def location_entries():
    return LocationEntry(
        ['chr1', 'chr2', 'chr2'],
        [19, 20, 50])

@pytest.fixture
def genomic_locations(genome_context):
    return GenomicLocation.from_fields(
        genome_context,
        ['chr1', 'chr2', 'chr2'],
        [19, 20, 50],
        ['-', '+', '+'])


@pytest.fixture
def track_arrays(chrom_sizes):
    return {name: np.repeat(np.arange(20), 5)[:size]
            for name, size in chrom_sizes.items()}


@pytest.fixture
def track(track_arrays, genome_context):
    return GenomicArrayGlobal.from_dict(
        {name: GenomicRunLengthArray.from_array(a)
         for name, a in track_arrays.items()}, genome_context)


@pytest.fixture
def track_stream(track_arrays, genome_context):
    return GenomicArrayNode.from_dict(
        {name: GenomicRunLengthArray.from_array(a)
         for name, a in track_arrays.items()}, genome_context)


@pytest.fixture
def genomic_intervals_big(intervals_big, genome_context_big):
    return GenomicIntervals.from_intervals(intervals_big, genome_context_big, is_stranded=False)


@pytest.fixture
def genomic_intervals(stranded_intervals, genome_context):
    return GenomicIntervals.from_intervals(stranded_intervals, genome_context, is_stranded=True)


@pytest.fixture
def genomic_intervals2(genome_context):
    return GenomicIntervals.from_fields(genome_context,
                                        ['chr1', 'chr1', 'chr2'],
                                        [5, 10, 0],
                                        [15, 95, 50])
