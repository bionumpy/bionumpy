import pytest

import bionumpy.encoded_array
from bionumpy.encoded_array import as_encoded_array
from bionumpy.encodings import StrandEncoding
from bionumpy.io.motifs import Motif
from bionumpy.simulate.chipseq import simulate_chip_seq_fragments, simulate_read_fragments
from bionumpy.simulate.rnaseq import get_transcript_copies, fragment_transcript_copies, sample_transcript_fragments, \
    get_rnaseq_reads, RNASeqSimulationSettings, simulate_rnaseq
import bionumpy as bnp
import numpy as np
from itertools import chain
from numpy.random import default_rng

rng = default_rng()


@pytest.fixture
def motif():
    return Motif("acgt",
                 np.array([[10, 15, 20, 10],
                           [12, 2, 10, 10],
                           [8, 13, 0, 10],
                           [5, 5, 5, 5]], dtype=int))


@pytest.fixture
def sequence():
    return bionumpy.encoded_array.as_encoded_array("acgtgcgtagctggctagctgcttagctgatggcttcgaa")


def test_simulate_chipseq(sequence, motif):
    fragments = simulate_chip_seq_fragments(sequence, motif, n_fragments=100, fragment_size=10)
    assert isinstance(fragments, bnp.datatypes.Interval)
    reads = simulate_read_fragments(fragments, 5)
    assert isinstance(reads, bnp.datatypes.Interval)


@pytest.fixture
def sequences():
    return bionumpy.encoded_array.as_encoded_array(["ACGT", "GCTA", "GTAAAT"], bnp.DNAEncoding)


@pytest.fixture
def sequence_counts():
    return [2, 1, 3]


@pytest.fixture
def rnaseq_simulation_settings(sequence_counts):
    return RNASeqSimulationSettings(fragment_size=3, read_length=2, transcript_counts=sequence_counts)


def test_get_transcript_copies(sequences, sequence_counts):
    print("SEQUENCES: %s" % repr(sequences[0].to_string()))
    truth = bionumpy.encoded_array.as_encoded_array(
        list(chain(*[[sequence.to_string()] * count for sequence, count in zip(sequences, sequence_counts)])),
        bnp.DNAEncoding)
    result = get_transcript_copies(sequences, sequence_counts)
    bnp.util.testing.assert_encoded_raggedarray_equal(truth, result)


def test_fragment_transcript_copies(sequences, fragment_size=2):
    truth = []
    for sequence in sequences:
        for i in range(0, len(sequence) - fragment_size + 1, fragment_size):
            truth.append(sequence[i:i + fragment_size])
    result = fragment_transcript_copies(sequences, fragment_size)
    bnp.util.testing.assert_encoded_raggedarray_equal(
        bionumpy.encoded_array.as_encoded_array(truth, sequences.encoding), bionumpy.encoded_array.as_encoded_array(result, sequences.encoding))


def test_sample_transcript_fragments(sequences, sampling_rate=0.9):
    np.random.seed(seed=123)
    mask = np.random.choice(a=[True, False], size=len(sequences), p=[sampling_rate, 1 - sampling_rate])
    truth = sequences[mask]
    result = sample_transcript_fragments(sequences, sampling_rate=0.9)
    bnp.util.testing.assert_encoded_raggedarray_equal(bionumpy.encoded_array.as_encoded_array(truth), bionumpy.encoded_array.as_encoded_array(result))


def test_get_rnaseq_reads(sequences, read_length=3):
    pos_strands = as_encoded_array("+" * len(sequences), StrandEncoding)
    pos_truth = as_encoded_array(["ACG", "GCT", "GTA"], sequences.encoding)
    neg_strands = as_encoded_array("-" * len(sequences), StrandEncoding)
    neg_truth = as_encoded_array(["ACG", "TAG", "ATT"], sequences.encoding)
    pos_result = get_rnaseq_reads(sequences, read_length, strands=pos_strands)
    neg_result = get_rnaseq_reads(sequences, read_length, strands=neg_strands)
    bnp.util.testing.assert_encoded_raggedarray_equal(pos_truth, pos_result)
    bnp.util.testing.assert_encoded_raggedarray_equal(neg_truth, neg_result)


def test_simualte_rnaseq(sequences, rnaseq_simulation_settings):
    result = simulate_rnaseq(sequences, rnaseq_simulation_settings)
    assert np.all([len(item.sequence) for item in result] == [rnaseq_simulation_settings.read_length] * len(result))

