from bionumpy.minimizers import get_minimizers
import numpy as np
import pytest

@pytest.fixture
def sequence():
    return np.array([0, 3, 1, 2, 2, 1, 0])


def test_minimizers(sequence):
    minimizers = get_minimizers(sequence, 4, 2)
    #assert np.testing.assert_equal(minizers
