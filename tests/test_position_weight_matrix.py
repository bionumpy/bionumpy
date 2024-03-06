import pytest
import numpy as np
from npstructures.testing import assert_raggedarray_equal
from numpy.testing import assert_array_equal
import bionumpy as bnp
from bionumpy.io.jaspar import read_jaspar_matrix
from bionumpy.sequence.position_weight_matrix import PositionWeightMatrix, PWM, get_motif_scores, get_motif_scores_old
from bionumpy import EncodedArray
from bionumpy.io.motifs import read_motif


@pytest.fixture
def neutral_ppm_dict():
    return {"A": [0.25, 0.25],
            "C": [0.25, 0.25],
            "G": [0.25, 0.25],
            "T": [0.25, 0.25]}


@pytest.fixture
def a_ppm_dict():
    return {"A": [1, 1],
            "C": [0, 0],
            "G": [0, 0],
            "T": [0, 0]}


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
    np.testing.assert_allclose(np.exp(log_prob), 0.4 * 0.25)


def test_sequence(sequence, matrix):
    pwm = PWM(matrix, "ACGT")
    log_prob = PositionWeightMatrix(pwm).rolling_window(sequence)
    np.testing.assert_allclose(np.exp(log_prob), [0.4 * 0.25, 0.025, 0.4 * 0.25])


def test_integration(data_path):
    # Read the alphabet and counts from jaspar file
    pwm = read_jaspar_matrix(data_path /"MA0080.1.jaspar")

    # Convert counts to position weight matrix
    # pwm = PWM.from_dict(pwm)

    # Make an array-class for the alphabet
    # encoding = AlphabetEncoding(alphabet)

    # Get the motif score function
    # pwm = PWM(pwm, alphabet)
    motif_score = PositionWeightMatrix(pwm)

    # Get reads
    entries = bnp.open(data_path / "reads.fq").read()

    # Calculate the motif score for each valid window
    scores = motif_score.rolling_window(entries.sequence)


def test_read_csv_motif(data_path):
    pwm = read_motif(data_path / "pwm.csv")
    pwm_jaspar = read_motif(data_path / "pwm.jaspar")
    assert str(pwm) == str(pwm_jaspar)


def test_pwm(window, matrix):
    pwm = PWM(matrix, "ACGT")
    # window = EncodedArray(window, AlphabetEncoding("ACGT"))
    log_prob = pwm.calculate_score(window)
    np.testing.assert_allclose(np.exp(log_prob), 0.4 * 0.25)


def test_encoded_ragged_array(sequences, matrix):
    pwm = PWM(matrix, "ACGT")
    get_motif_scores(sequences, pwm)


def test_encoded_ragged_array_fast(sequences, matrix):
    pwm = PWM(matrix, "ACGT")
    s = get_motif_scores_old(sequences, pwm)
    assert_raggedarray_equal(s,
                             get_motif_scores(sequences, pwm))


def test_sanity_motifs(sequence, neutral_ppm_dict):
    pwm = PWM.from_dict(neutral_ppm_dict)
    scores = pwm.calculate_scores(sequence)
    assert np.all(scores == 0)


def test_a_motifs(a_ppm_dict):
    pwm = PWM.from_dict(a_ppm_dict)
    scores = pwm.calculate_scores("AAC")
    assert_array_equal(scores, [np.log(4 ** 2), -np.inf, -np.inf])


@pytest.mark.skip("Failing because under development?")
def test_from_dict(window, matrix):
    dictionary = dict(zip("ACGT", matrix))
    # window = EncodedArray(window, AlphabetEncoding("ACGT"))
    pwm = PWM.from_dict(dictionary)
    log_prob = pwm.calculate_score(window)
    np.testing.assert_allclose(np.exp(log_prob), 0.4 * 0.25)
