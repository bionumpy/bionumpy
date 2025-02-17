import pytest

import bionumpy.encoded_array
from bionumpy.encodings.kmer_encodings import KmerEncoding
from bionumpy.sequence.kmers import KmerEncoder, _get_dna_kmers
import numpy as np
from npstructures import RaggedShape
from bionumpy import DNAEncoding, EncodedRaggedArray, EncodedArray, as_encoded_array
import bionumpy as bnp


def test_kmer_encoding():
    encoding = KmerEncoder(3, DNAEncoding)
    kmers = encoding.sample_domain(100)
    encoded = encoding(kmers)
    decoded = encoding.inverse(encoded)
    np.testing.assert_equal(np.asarray(kmers), np.asarray(decoded))


@pytest.mark.parametrize("sequences", [
    "ACTG",
    "ACACATCGACGAgactagct",
    "AacACtggatcggacTTATCTGACG",
    "G"
])

def test_get_dna_kmers_equals_get_kmers(sequences):
    kmers = _get_dna_kmers(as_encoded_array("cgtt", DNAEncoding), 3)
    encoding = KmerEncoder(3, DNAEncoding)
    np.testing.assert_array_equal(kmers, encoding.rolling_window("cgtt"))


def test_rolling_hash():
    lengths = np.arange(3, 10)
    encoding = KmerEncoder(3, DNAEncoding)
    kmers = EncodedArray(np.arange(lengths.sum()) % 4, DNAEncoding)
    ragged = EncodedRaggedArray(kmers, lengths)
    encoded = encoding.rolling_window(ragged)
    encoded.ravel()
    assert encoded._shape == RaggedShape(lengths-3+1)


def test_get_kmers():
    sequence = bionumpy.encoded_array.as_encoded_array(
        ["ACTG", "CAAAAA", "TTT"], bnp.DNAEncoding
    )
    correct = [
        ["ACT", "CTG"],
        ["CAA", "AAA", "AAA", "AAA"],
        ["TTT"]
    ]
    kmers = bnp.sequence.get_kmers(sequence, 3)
    decoded_kmers = [
        [str(k) for k in read_kmers] for read_kmers in kmers
    ]

    print(decoded_kmers == correct)

@pytest.mark.parametrize("encoding", [bnp.DNAEncoding, bnp.AminoAcidEncoding])
def test_get_kmers_one(encoding):
    sequence = bionumpy.encoded_array.as_encoded_array(["ACTG"], encoding)
    kmers = bnp.sequence.get_kmers(sequence, 1)
    assert len(kmers[0]) == 4, kmers[0]


@pytest.mark.skip("Not correct")
def test_kmer_encoding():
    encoding = bnp.encodings.kmer_encodings.KmerEncoding(
        bnp.encodings.DNAEncoding, k=5
    )
    assert encoding.to_string(1) == "AAAAC"
    assert encoding.to_string(2) == "AAAAG"
    assert encoding.to_string(4) == "AAACA"


def test_kmer_encoding_repr_and_string():
    sequences = bionumpy.encoded_array.as_encoded_array(["ACTG", "AAA", "TTGGC"], bnp.DNAEncoding)
    kmers = bnp.sequence.get_kmers(sequences, 3)
    print(str(kmers[0]))
    print(repr(kmers[0]))
    print(repr(kmers))
    print(str(kmers))


def test_get_kmer_encoding_labels():
    encoding = KmerEncoding(bnp.DNAEncoding, 3)
    labels = encoding.get_labels()
    assert labels[0:5] == [
        "AAA",
        "CAA",
        "GAA",
        "TAA",
        "ACA"
    ]


def test_count_kmers():
    sequences = bionumpy.encoded_array.as_encoded_array(["ACTG", "AAA", "TTGGC"], bnp.DNAEncoding)
    kmers = bnp.sequence.get_kmers(sequences, 3)
    counts = bnp.count_encoded(kmers, axis=None)
    assert counts["ACT"] == 1
    assert counts["GGG"] == 0

