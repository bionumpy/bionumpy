import pytest
import numpy as np

import bionumpy as bnp
from bionumpy.io.jaspar import read_jaspar_matrix
from bionumpy.sequence.position_weight_matrix import PositionWeightMatrix, pwm_from_counts
from bionumpy.encodings.alphabet_encoding import AlphabetEncoding


@pytest.fixture
def matrix():
    with np.errstate(divide='ignore'):
        m = np.log([[0.4, 0.25],
                    [0.1, 0.25],
                    [0.4, 0.25],
                    [0.1, 0.25]])
    return m


@pytest.fixture
def window():
    return np.array([0, 1])


@pytest.fixture
def sequence():
    return np.array([0, 1, 2, 3])


def test_window(window, matrix):
    log_prob = PositionWeightMatrix(matrix)(window)
    np.testing.assert_allclose(np.exp(log_prob), 0.4*0.25)


def test_sequence(sequence, matrix):
    log_prob = PositionWeightMatrix(matrix).rolling_window(sequence)
    np.testing.assert_allclose(np.exp(log_prob), [0.4*0.25, 0.025, 0.4*0.25])


def test_integration():
    # Read the alphabet and counts from jaspar file
    alphabet, matrix = read_jaspar_matrix("example_data/MA0080.1.jaspar")
    
    # Convert counts to position weight matrix
    pwm = pwm_from_counts(matrix)
    
    # Make an array-class for the alphabet
    encoding = AlphabetEncoding(alphabet)
    
    # Get the motif score function
    motif_score = PositionWeightMatrix(pwm, encoding)
    
    #Get reads
    entries = bnp.open("example_data/reads.fq").read()
    
    # Calculate the motif score for each valid window
    scores = motif_score.rolling_window(entries.sequence)
    
