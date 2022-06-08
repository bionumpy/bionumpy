from bionumpy.kmers import KmerEncoding
import numpy as np
from npstructures import RaggedArray, RaggedShape


def test_kmer_encoding():
    encoding = KmerEncoding(3, alphabet_size=4)
    kmers = encoding.sample_domain(100)
    encoded = encoding(kmers)
    # assert np.testing.assert_ encoding.in_range(encoded)
    decoded = encoding.inverse(encoded)
    np.testing.assert_equal(kmers, decoded)


def test_rolling_hash():
    lengths = np.arange(3, 10)
    encoding = KmerEncoding(3, alphabet_size=4)
    kmers = np.arange(lengths.sum()) % 4
    ragged = RaggedArray(kmers, lengths)
    encoded = encoding.rolling_window(ragged)
    assert encoded.shape == RaggedShape(lengths-3+1)
