from bionumpy import as_encoded_array, DNAEncoding
from bionumpy.sequence.indexing import KmerIndex, KmerLookup
import pytest
import numpy as np


@pytest.fixture
def sequences():
    return as_encoded_array(["ACGTAA", 'GCTAAA'], DNAEncoding)


def test_kmer_index(sequences):
    index = KmerIndex.create_index(sequences, k=3)
    assert index.get_indices("ACG") == [0]
    assert index.get_indices("AAA") == [1]
    np.testing.assert_equal(index.get_indices("TAA"), [0, 1])
    assert index.get_indices("GAA") == []


def test_multiple_occurrences(sequences):
    index = KmerIndex.create_index(sequences, k=2)
    np.testing.assert_equal(index.get_indices("AA"), [0, 1])


def test_kmer_lookup(sequences):
    lookup = KmerLookup.create_lookup(sequences, k=3)
    sequences_with_kmer = lookup.get_sequences(kmer="CGT")
    assert sequences_with_kmer == ["ACGTAA"]
