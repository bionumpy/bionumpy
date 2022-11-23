import pytest
from bionumpy.io.motifs import Motif, read_motif
from bionumpy.simulate.chipseq import simulate_chip_seq_fragments, simulate_read_fragments
from bionumpy.simulate.rnaseq import get_transcript_copies, fragment_transcript_copies
from bionumpy.sequence import get_kmers
import bionumpy as bnp
import numpy as np
from itertools import chain


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

@pytest.fixture
def sequences():
    return bnp.as_encoded_array(["ACGT", "GCTA", "GTAAAT"], bnp.DNAEncoding)

@pytest.fixture
def sequence_counts():
    return [2,1,3]

def test_get_transcript_copies(sequences, sequence_counts):
    truth = bnp.as_encoded_array(list(chain(*[[sequence.to_string()]*count for sequence,count in zip(sequences, sequence_counts)])), bnp.DNAEncoding)
    result = get_transcript_copies(sequences, sequence_counts)
    bnp.testing.assert_encoded_raggedarray_equal(truth,result)

def test_fragment_transcript_copies(sequences, fragment_size=2):
    truth = []
    for sequence in sequences:
        for i in range(0, len(sequence) - fragment_size + 1, fragment_size):
            truth.append(sequence[i:i + fragment_size])
    result = fragment_transcript_copies(sequences, fragment_size)
    bnp.testing.assert_encoded_raggedarray_equal(truth,result)

def test_sample_transcript_fragments(sequences, sampling_rate=0.9):
    pass

