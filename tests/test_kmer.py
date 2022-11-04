from bionumpy.kmers import KmerEncoding, fast_hash
import numpy as np
from npstructures import RaggedShape
from bionumpy import DNAEncoding, EncodedRaggedArray, EncodedArray, as_encoded_array


def test_kmer_encoding():
    encoding = KmerEncoding(3, DNAEncoding)
    kmers = encoding.sample_domain(100)
    encoded = encoding(kmers)
    decoded = encoding.inverse(encoded)
    np.testing.assert_equal(np.asarray(kmers), np.asarray(decoded))


def test_fast_hash():
    kmers = fast_hash(as_encoded_array("cgtt", DNAEncoding), 3)
    encoding = KmerEncoding(3, DNAEncoding)
    np.testing.assert_array_equal(kmers, encoding.rolling_window("cgtt"))


def test_rolling_hash():
    lengths = np.arange(3, 10)
    encoding = KmerEncoding(3, DNAEncoding)
    kmers = EncodedArray(np.arange(lengths.sum()) % 4, DNAEncoding)
    ragged = EncodedRaggedArray(kmers, lengths)
    encoded = encoding.rolling_window(ragged)
    assert encoded.shape == RaggedShape(lengths-3+1)
