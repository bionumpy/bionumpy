import pytest
import bionumpy as bnp
import numpy as np

import bionumpy.encoded_array
from bionumpy import as_encoded_array
from bionumpy.sequence.indexing import WildCardIndex, WildCardLookup
from bionumpy.util.testing import assert_encoded_raggedarray_equal


@pytest.fixture
def sequences():
    return bionumpy.encoded_array.as_encoded_array(["ACGTACG", "CGTAGT"], bnp.DNAEncoding)


def test_wildcard_index(sequences):
    index = WildCardIndex.create_index(sequences)
    np.testing.assert_equal(index.get_indices(pattern="A.G"), [0])
    np.testing.assert_equal(index.get_indices("CGT"), [0, 1])
    np.testing.assert_equal(index.get_indices("ACG"), [0])
    np.testing.assert_equal(index.get_indices("ACGCG"), [])
    np.testing.assert_equal(index.get_indices("A..TA"), [0])


def test_wildcard_lookup(sequences):
    lookup = WildCardLookup.create_lookup(sequences)
    assert lookup.get_sequences("TA.") == ["ACGTACG", "CGTAGT"]
    assert len(lookup.get_sequences("GC.")) == 0


