from bionumpy.kmers import KmerEncoding, safe_hash
import numpy as np
from npstructures import RaggedArray, RaggedShape
from bionumpy.encodings import ACTGEncoding


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



def test_safe_hash():
    sequence = np.array([1, 2, 3, 1, 2, 3, 1, 1, 1, 1], dtype=np.uint8)

    k = 5
    hashes = safe_hash(sequence, k)
    print(hashes)

    hasher = KmerEncoding(k, None, 4)
    hashes = hasher.rolling_window(sequence)
    #print(hasher(sequence))
    print(type(hashes))
    print(hashes)
