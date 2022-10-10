from bionumpy.kmers import KmerEncoding
import numpy as np
from npstructures import RaggedArray, RaggedShape
from bionumpy.encodings import ACTGEncoding
from bionumpy import DNAArray


def test_kmer_encoding():
    encoding = KmerEncoding(3, DNAArray)
    kmers = encoding.sample_domain(100)
    encoded = encoding(kmers)
    # assert np.testing.assert_ encoding.in_range(encoded)
    decoded = encoding.inverse(encoded)
    print(kmers, decoded)
    np.testing.assert_equal(np.asarray(kmers), np.asarray(decoded))


def test_rolling_hash():
    lengths = np.arange(3, 10)
    encoding = KmerEncoding(3, DNAArray)
    kmers = np.arange(lengths.sum()) % 4
    ragged = RaggedArray(kmers, lengths)
    encoded = encoding.rolling_window(ragged)
    assert encoded.shape == RaggedShape(lengths-3+1)

