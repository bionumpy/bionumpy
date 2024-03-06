import pytest
import bionumpy as bnp
from bionumpy.encodings.string_encodings import StringEncoding, AsciiHashTable
from bionumpy.encodings.exceptions import EncodingError
from numpy.testing import assert_array_equal
import bionumpy as bnp


@pytest.fixture
def encoded_ragged_array():
    return bnp.as_encoded_array(["chr1",
                                 "chr2",
                                 "chr13"])

@pytest.fixture
def problems():
    return bnp.as_encoded_array([
        'chr4_GL000257v2_alt',
        'chr6_GL000256v2_alt'])


@pytest.fixture
def string_encoding(encoded_ragged_array):
    return StringEncoding(encoded_ragged_array, modulo=103)

@pytest.fixture
def chrom_names(data_path):
    return bnp.open(data_path / "hg38.chrom.sizes").read().name


def test_problems(problems):
    AsciiHashTable.from_sequences(problems, modulo=103)


def test_chromosome_encoding(chrom_names):
    from collections import Counter
    print(Counter(chrom_names.tolist()).most_common(10))
    encoding = AsciiHashTable.from_sequences(chrom_names, modulo=103)


def test_string_encoding_labels(encoded_ragged_array):
    encoding = StringEncoding(encoded_ragged_array, modulo=103)
    assert encoding.get_labels() == ["chr1", "chr2", "chr13"]


def test_string_encoding_encode(encoded_ragged_array):
    encoding = StringEncoding(encoded_ragged_array, modulo=103)
    encoded = encoding.encode(encoded_ragged_array)
    assert_array_equal(encoded.raw(), [0, 1, 2])


def test_string_encoding_encode_invalid(string_encoding, problems):
    with pytest.raises(EncodingError):
        string_encoding.encode(problems)
        # encoding = StringEncoding(encoded_ragged_array, modulo=103)
#    encoded = encoding.encode(encoded_ragged_array)
#     assert_array_equal(encoded.raw(), [0, 1, 2])
