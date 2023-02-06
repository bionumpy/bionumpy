from bionumpy.genomic_data import GenomicTrack
from bionumpy.streams import NpDataclassStream
from bionumpy.datatypes import BedGraph
import numpy as np
from numpy.testing import assert_equal
import pytest


@pytest.fixture
def chrom_sizes():
    return {"chr1": 100, "chr2": 50}


@pytest.fixture
def pileup():
    return BedGraph.from_entry_tuples([
        ('chr1', 0, 10, 1),
        ('chr1', 10, 90, 2),
        ('chr1', 90, 100, 3),
        ('chr2', 0, 50, 4)])


@pytest.fixture
def track(pileup, chrom_sizes):
    return GenomicTrack.from_bedgraph(pileup, chrom_sizes)


@pytest.fixture
def stream_track(pileup, chrom_sizes):
    return GenomicTrack.from_bedgraph(NpDataclassStream([pileup]), chrom_sizes)


@pytest.fixture
def array(track):
    return np.concatenate([rle for rle in track.to_dict().values()])


def test_extract_chromosome(track):
    assert_equal(track['chr2'].to_array(),
                 np.full(50, 4))


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
