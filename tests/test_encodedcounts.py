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


def test_ufuncs(encoded_counts):
    sqrt = np.sqrt(encoded_counts)
    assert sqrt == EncodedCounts(list('ACGT'), np.sqrt([4, 3, 1, 5]))
    assert encoded_counts + 1 == EncodedCounts(list('ACGT'), np.array([5, 4, 2, 6]))
    assert encoded_counts + encoded_counts == EncodedCounts(list('ACGT'), np.array([8, 6, 2, 10]))

