from npstructures import RaggedArray
from npstructures.testing import assert_raggedarray_equal
from bionumpy.genomic_data import GenomicArray, GenomicIntervals
from bionumpy.genomic_data.genome_context import GenomeContext
from bionumpy.streams import NpDataclassStream
from bionumpy.datatypes import BedGraph
import numpy as np
from numpy.testing import assert_equal
import pytest


@pytest.fixture
def chrom_sizes():
    return {"chr1": 100, "chr2": 50}


@pytest.fixture
def genome_context(chrom_sizes):
    return GenomeContext.from_dict(chrom_sizes)


@pytest.fixture
def pileup():
    return BedGraph.from_entry_tuples([
        ('chr1', 0, 10, 1),
        ('chr1', 10, 90, 2),
        ('chr1', 90, 100, 3),
        ('chr2', 0, 50, 4)])


@pytest.fixture
def genomic_intervals(genome_context):
    return GenomicIntervals.from_fields(genome_context,
                                        ['chr1', 'chr1', 'chr2'],
                                        [5, 10, 0],
                                        [15, 95, 50],
                                        ['-', '+', '-'])


@pytest.fixture
def interval_pileup():
    return RaggedArray(
        [np.concatenate([np.full(5, 2), np.full(5, 1)]),
         np.concatenate([np.full(80, 2), np.full(5, 3)]),
         np.full(50, 4)])


@pytest.fixture
def track(pileup, genome_context) -> GenomicArray:
    return GenomicArray.from_bedgraph(pileup, genome_context)


@pytest.fixture
def stream_track(pileup, genome_context):
    return GenomicArray.from_bedgraph(NpDataclassStream([pileup]), genome_context)


@pytest.fixture
def array(track):
    return np.concatenate([rle for rle in track.to_dict().values()])


def test_extract_chromosome(track):
    assert_equal(track['chr2'].to_array(),
                 np.full(50, 4))


def test_extract_intervals(track, genomic_intervals, interval_pileup):
    result = track[genomic_intervals]
    assert_raggedarray_equal(result.to_array(),
                             interval_pileup)


def test_histogram(track, array):
    hist = np.histogram(track)
    true_hist = np.histogram(array)
    for h, t in zip(hist, true_hist):
        np.testing.assert_equal(h, t)


def test_stream_histogram(stream_track, array):
    hist = np.histogram(stream_track, bins=3, range=(0, 5)).compute()
    true_hist = np.histogram(array, bins=3, range=(0, 5))
    for h, t in zip(hist, true_hist):
        np.testing.assert_equal(h, t)


