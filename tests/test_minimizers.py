import bionumpy.encoded_array
from bionumpy.encodings.kmer_encodings import KmerEncoding
from bionumpy.sequence.minimizers import Minimizers, get_minimizers
from bionumpy.encoded_array import EncodedArray, EncodedRaggedArray
from bionumpy.sequence.kmers import KmerEncoder
from bionumpy.encodings.alphabet_encoding import DNAEncoding
from npstructures import RaggedArray
import numpy as np
import pytest
import bionumpy as bnp
from bionumpy.util.testing import assert_encoded_raggedarray_equal, assert_encoded_array_equal

@pytest.fixture
def sequence():
    return EncodedArray(np.array([0, 3, 1, 2, 2, 1, 0]), DNAEncoding)


@pytest.fixture
def sequences():
    r = RaggedArray([[0, 3, 1, 2, 2, 1, 0],
                     [0, 3, 1, 2, 2, 1],
                     [0, 3, 1, 2, 2],
                     [0, 3, 1, 2]])
    r = EncodedRaggedArray(
        EncodedArray(r.ravel(), DNAEncoding), r._shape)
    return r


@pytest.fixture
def window():
    return EncodedArray(np.array([0, 3, 1, 2]), DNAEncoding)


@pytest.fixture
def encoding():
    return Minimizers(3, KmerEncoder(2, DNAEncoding))


@pytest.fixture()
def kmer_encoding():
    return KmerEncoding(DNAEncoding, 2)


def test_window(window, encoding, kmer_encoding):
    minimizer = encoding(window)
    true = EncodedArray([7], kmer_encoding)
    np.testing.assert_equal(minimizer, true)


def test_minimizers(sequence, encoding, kmer_encoding):
    minimizers = get_minimizers(sequence, 2, 4)
    true = EncodedArray([7, 7, 6, 1], kmer_encoding)
    assert_encoded_array_equal(minimizers, true)


def test_full_roll(sequences, encoding, kmer_encoding):
    k = 2
    w = 4
    minimizers = get_minimizers(sequences, k, w)
    ra = RaggedArray([[7, 7, 6, 1], [7, 7, 6], [7, 7], [7]])
    true = EncodedRaggedArray(EncodedArray(ra.ravel(), kmer_encoding), ra._shape)
    assert_encoded_raggedarray_equal(true, minimizers)


def test_get_minimizer_string_to_string():
    sequences = bionumpy.encoded_array.as_encoded_array([
        "CCCAAACCCC",
        "TTTTCCCTTT"
    ], DNAEncoding)

    k = 3
    w = 10

    minimizers = get_minimizers(sequences, k, w)
    string_minimizers = [
        [str(m) for m in sequence_minimizers]
        for sequence_minimizers in minimizers]

    correct = [["AAA"], ["CCC"]]
    assert string_minimizers == correct

