import pytest
import numpy as np

from bionumpy.sequence.bloom_filter import BloomFilter, hash_function


@pytest.fixture
def kmer_hashes():
    return np.arange(0, 64, 5)


@pytest.fixture
def hash_functions():
    return [hash_function(offset, 17) for offset in (13, 51)]


def test_bloom_filter_precence(kmer_hashes, hash_functions):
    f = BloomFilter.from_hash_functions_and_seqeuences(hash_functions, kmer_hashes, 17)
    assert np.all(f[kmer_hashes])
