import pytest
import numpy as np

import bionumpy as bnp
from bionumpy.io.jaspar import read_jaspar_matrix
from bionumpy.sequence.position_weight_matrix import PositionWeightMatrix, _pwm_from_counts, PWM, get_motif_scores
from bionumpy.encodings.alphabet_encoding import AlphabetEncoding
from bionumpy import EncodedArray


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
    return EncodedArray(np.array([0, 1]), bnp.DNAEncoding)


@pytest.fixture
def sequence():
    return EncodedArray(np.array([0, 1, 2, 3]), bnp.DNAEncoding)


@pytest.fixture
def sequences():
    return ["ACGT",
            "GCT"]


def test_window(window, matrix):
    pwm = PWM(matrix, "ACGT")
    log_prob = PositionWeightMatrix(pwm)(window)
    np.testing.assert_allclose(np.exp(log_prob), 0.4*0.25)


def test_sequence(sequence, matrix):
    pwm = PWM(matrix, "ACGT")
    log_prob = PositionWeightMatrix(pwm).rolling_window(sequence)
    np.testing.assert_allclose(np.exp(log_prob), [0.4*0.25, 0.025, 0.4*0.25])


def test_integration():
    # Read the alphabet and counts from jaspar file
    alphabet, matrix = read_jaspar_matrix("example_data/MA0080.1.jaspar")
    
    # Convert counts to position weight matrix
    pwm = _pwm_from_counts(matrix)
    
    # Make an array-class for the alphabet
    encoding = AlphabetEncoding(alphabet)
    
    # Get the motif score function
    pwm = PWM(pwm, alphabet)
    motif_score = PositionWeightMatrix(pwm)
    
    #Get reads
    entries = bnp.open("example_data/reads.fq").read()
    
    # Calculate the motif score for each valid window
    scores = motif_score.rolling_window(entries.sequence)
    

def test_pwm(window, matrix):
    pwm = PWM(matrix, "ACGT")
    #window = EncodedArray(window, AlphabetEncoding("ACGT"))
    log_prob = pwm.calculate_score(window)
    np.testing.assert_allclose(np.exp(log_prob), 0.4*0.25)


def test_encoded_ragged_array(sequences, matrix):
    pwm = PWM(matrix, "ACGT")
    get_motif_scores(sequences, pwm)


@pytest.mark.skip("Failing because under development?")
def test_from_dict(window, matrix):
    dictionary = dict(zip("ACGT", matrix))
    # window = EncodedArray(window, AlphabetEncoding("ACGT"))
    pwm = PWM.from_dict(dictionary)
    log_prob = pwm.calculate_score(window)
    np.testing.assert_allclose(np.exp(log_prob), 0.4*0.25)
    
