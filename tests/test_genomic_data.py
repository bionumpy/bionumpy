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

from .genomic_fixtures import *


def test_get_location(genomic_intervals, genomic_locations):
    result = genomic_intervals.get_location('start')
    assert_equal(result.position, genomic_locations.position)


def test_map_locations(genomic_intervals2):
    locations = LocationEntry(['chr1',
                               'chr1',
                               'chr2'], [5, 10, 0])
    result = genomic_intervals2.map_locations(locations)
    true_mapped = LocationEntry(['0', '0', '1', '2'],
                                [0, 5, 0, 0])
    assert_bnpdataclass_equal(result, true_mapped)


def test_get_windows(genomic_locations, windows):
    result = genomic_locations.get_windows(flank=5)
    assert_equal(result.start, windows.start)
    assert_equal(result.stop, windows.stop)


def test_extract_windows(windows, track_stream, track_arrays):
    tmp = track_stream[windows].mean(axis=0)
    result = compute(tmp)
    result = result.to_array()
    true = np.mean([track_arrays[w.chromosome.to_string()][w.start:w.stop] if w.strand == '+' else track_arrays[w.chromosome.to_string()][w.stop-1:w.start-1:-1] for w in windows], axis=0)
    assert_raggedarray_equal(result, true)


def test_extract_windows_2(windows, track_stream, track_arrays):
    tmp = track_stream[windows]
    result = compute(tmp)
    result = result.to_array()
    true = np.array([track_arrays[w.chromosome.to_string()][w.start:w.stop] if w.strand == '+' else track_arrays[w.chromosome.to_string()][w.stop-1:w.start-1:-1] for w in windows])
    assert_equal(result, true)


def test_extract_windows_and_subset(windows, track_stream, track_arrays):
    tmp = track_stream[windows]
    tmp = tmp[:, ::2]
    result = compute(tmp)
    result = result[np.argsort(windows.stop-windows.start)]

    result = result.to_array()
    true = np.array([track_arrays[w.chromosome.to_string()][w.start:w.stop][::2]
                    if w.strand == '+' else track_arrays[w.chromosome.to_string()][w.stop-1:w.start-1:-1][::2] for w in windows])
    assert_equal(result, true)


def test_genome_context(genome_context_big):
    locations = LocationEntry(
        ['chr1', 'chr10_alt', 'chr3'],
        [0, 1, 2])
    genome_context = genome_context_big
    filled = genome_context.iter_chromosomes(locations, LocationEntry)
    true = [LocationEntry.from_entry_tuples([('chr1', 0)]),
            LocationEntry.empty(),
            LocationEntry.from_entry_tuples([('chr3', 2)])]
    for f, t in zip(filled, true):
        assert_bnpdataclass_equal(f, t)


def test_global_offset(genome_context_big, genomic_intervals_big):
    genomic_intervals_big.get_pileup()


def test_empty_intervals_pileup(genome_context_big, genomic_intervals_big):
    genomic_intervals_big.as_stream().get_pileup()


def test_invalid_context(genome_context_big):
    genome_context = genome_context_big
    locations = LocationEntry(
        ['chr3', 'chr10_alt', 'chr1'],
        [0, 1, 2])
    with pytest.raises(Exception):
        list(genome_context.iter_chromosomes(locations, LocationEntry))


def test_invalid_context2(genome_context_big):
    genome_context = genome_context_big
    locations = LocationEntry(
        ['chr1', 'chr11_alt', 'chr3'],
        [0, 1, 2])
    with pytest.raises(GenomeError):
        list(genome_context.iter_chromosomes(locations, LocationEntry))


