import pytest
import bionumpy as bnp
import numpy as np


@pytest.fixture
def array_1():
    return bnp.as_encoded_array("ACGTTGCA", target_encoding=bnp.DNAEncoding)

@pytest.fixture
def array_2():
    return bnp.as_encoded_array("TTTTGGGG", target_encoding=bnp.DNAEncoding)


def test_lex_sort(array_1):
    np.testing.assert_array_equal(np.lexsort((array_1,)), [0, 7, 1, 6, 2, 5, 3, 4])


def test_lex_sort_2(array_1, array_2):
    np.testing.assert_array_equal(np.lexsort((array_1, array_2)), 
                                  [7, 6, 5, 4, 0, 1, 2, 3])
