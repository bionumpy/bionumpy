import pytest
import bionumpy as bnp
import numpy as np
from bionumpy.wildcard_index import WildCardIndex


@pytest.fixture
def sequences():
    return bnp.as_encoded_array(["ACGTACG", "CGTAGT"], bnp.DNAEncoding)


def test_wildcard_index(sequences):
    index = WildCardIndex.create_index(sequences)
    pattern_indices = index.get_pattern_indices(pattern="A.G")
    np.testing.assert_equal(pattern_indices, [0])

