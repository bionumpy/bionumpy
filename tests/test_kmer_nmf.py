import pytest
import numpy as np
import bionumpy as bnp
from bionumpy import as_encoded_array
#from scripts.kmer_nmf import estimate_reference_proportions


@pytest.fixture
def reference_sequences():
    reference_sequences = as_encoded_array(["AAAAAGGGGGCCCCCTTTTTT", "ACGTACGTACGTACGT"], bnp.DNAEncoding)
    return reference_sequences


@pytest.fixture
def query_sequences():
    return as_encoded_array(["ACG", "ACG", "ACG"], bnp.DNAEncoding)

@pytest.mark.skip
def test_estimate_reference_proportions(reference_sequences, query_sequences):
    proportions = estimate_reference_proportions(reference_sequences, query_sequences)
    np.testing.assert_array_equal(proportions, [0, 1])
