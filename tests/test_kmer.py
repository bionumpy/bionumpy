from bionumpy.kmers import KmerEncoding
import numpy as np


def test_kmer_encoding():
    encoding = KmerEncoding(3, 4)
    kmers = encoding.sample_domain(100)
    encoded = encoding.from_bytes(kmers)
    assert encoding.in_range(encoded)
    decoded = encoding.to_bytes(encoded)
    np.testing.assert_equal(kmers, decoded)


def test_hash_sequences():
    pass
