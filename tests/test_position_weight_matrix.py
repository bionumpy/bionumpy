import pytest
import numpy as np

from bionumpy.position_weight_matrix import PositionWeightMatrix


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
