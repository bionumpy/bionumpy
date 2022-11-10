import pytest
from bionumpy.io.motifs import Motif, read_motif
from bionumpy.simulate.chipseq import simulate_chip_seq_fragments, simulate_read_fragments
import bionumpy as bnp
import numpy as np


@pytest.fixture
def motif():
    return Motif("acgt",
                 np.array([[10, 15, 20, 10],
                           [12, 2, 10, 10],
                           [8, 13, 0, 10],
                           [5, 5, 5, 5]], dtype=int))


@pytest.fixture
def sequence():
    return bnp.as_encoded_array("acgtgcgtagctggctagctgcttagctgatggcttcgaa")


def test_simulate_chipseq(sequence, motif):
    fragments = simulate_chip_seq_fragments(sequence, motif, n_fragments=100, fragment_size=10)
    assert isinstance(fragments, bnp.datatypes.Interval)
    reads = simulate_read_fragments(fragments, 5)
    assert isinstance(reads, bnp.datatypes.Interval)
