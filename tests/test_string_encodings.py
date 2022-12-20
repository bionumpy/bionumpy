import pytest
from bionumpy.encodings.string_encodings import StringEncoding
from numpy.testing import assert_array_equal
import bionumpy as bnp


@pytest.fixture
def encoded_ragged_array():
    return bnp.as_encoded_array(["chr1",
                                 "chr2",
                                 "chr13"])


def test_string_encoding_labels(encoded_ragged_array):
    encoding = StringEncoding(encoded_ragged_array, modulo=103)
    assert encoding.get_labels() == ["chr1", "chr2", "chr13"]


def test_string_encoding_encode(encoded_ragged_array):
    encoding = StringEncoding(encoded_ragged_array, modulo=103)
    encoded = encoding.encode(encoded_ragged_array)
    assert_array_equal(encoded.raw(), [0, 1, 2])
