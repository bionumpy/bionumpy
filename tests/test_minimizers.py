from bionumpy.minimizers import Minimizers
from bionumpy.sequences import Sequences, Sequence
from bionumpy.kmers import KmerEncoding
from bionumpy.encodings import ACTGEncoding
import numpy as np
import pytest


@pytest.fixture
def sequence():
    return np.array([0, 3, 1, 2, 2, 1, 0])


@pytest.fixture
def sequences():
    return Sequences([[0, 3, 1, 2, 2, 1, 0],
                      [0, 3, 1, 2, 2, 1],
                      [0, 3, 1, 2, 2],
                      [0, 3, 1, 2]], encoding=ACTGEncoding)


@pytest.fixture
def window():
    s = Sequence.from_array(np.array([0, 3, 1, 2]))
    s.encoding = ACTGEncoding
    return s


@pytest.fixture
def encoding():
    return Minimizers(3, KmerEncoding(2, alphabet_size=4))


def test_window(window, encoding):
    minimizer = encoding(window)
    np.testing.assert_equal(minimizer, 7)


def test_minimizers(sequence, encoding):
    minimizers = encoding.rolling_window(sequence)
    np.testing.assert_equal(minimizers, [7, 7, 6, 1])


def test_full_roll(sequences, encoding):
    minimizers = encoding.rolling_window(sequences)
    true = [[7, 7, 6, 1], [7, 7, 6], [7, 7], [7]]
    for row, true_row in zip(minimizers, true):
        np.testing.assert_equal(row, true_row)
