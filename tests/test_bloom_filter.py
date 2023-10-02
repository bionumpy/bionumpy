import pytest
import numpy as np

from bionumpy.sequence.bloom_filter import BloomFilter, hash_function, InterleavedBloomFilter


@pytest.fixture
def kmer_hashes():
    return np.arange(0, 64, 5)


@pytest.fixture
def kmer_hashes_list():
    return np.arange(0, 64).reshape(8, 8)//np.arange(1, 9)[:, None]


@pytest.fixture
def hash_functions():
    return [hash_function(offset, 17) for offset in (13, 51)]

@pytest.mark.skip
def test_bloom_filter_precence(kmer_hashes, hash_functions):
    f = BloomFilter.from_hash_functions_and_seqeuences(hash_functions,
                                                       kmer_hashes, 17)
    assert np.all(f[kmer_hashes])


@pytest.mark.skip
def test_interleaved_bloom_filter_precence(kmer_hashes_list, hash_functions):
    f = InterleavedBloomFilter.from_hash_functions_and_seqeuences(hash_functions,
                                                                  kmer_hashes_list, 17)
    for i, seq in enumerate(kmer_hashes_list):
        assert np.all(f[seq, i])

@pytest.mark.skip
def test_needle_acceptance(sequences):
    minimizers = get_minimizers(sequences, k=31, w=41)

