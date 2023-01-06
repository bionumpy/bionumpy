from bionumpy.sequence.count_encoded import EncodedCounts
import numpy as np
import pytest


@pytest.fixture
def encoded_counts():
    return EncodedCounts(
        list('ACGT'),
        np.array([4, 3, 1, 5], dtype=int))


def test_most_common(encoded_counts):
    most_common = encoded_counts.most_common(2)
    assert most_common == EncodedCounts(
        ['T', 'A'], np.array([5, 4]))
