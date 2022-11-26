import numpy as np
import pytest
import bionumpy as bnp
from npstructures.testing import assert_raggedarray_equal
from npstructures import RaggedArray
import bionumpy.encoded_array
import bionumpy.encoded_array_functions
import bionumpy as bnp
from bionumpy.encodings import DigitEncoding
from bionumpy.encodings.base_encoding import NumericEncoding, OneToOneEncoding
from bionumpy.encoded_array_functions import as_encoded_array


#from bionumpy.encodings.base_encoding import OneToOneEncoding


def test_change_encoding_on_encoded_array():
    a = bionumpy.encoded_array_functions.as_encoded_array("ACTG", bnp.encodings.alphabet_encoding.ACTGEncoding)
    b = bionumpy.encoded_array_functions.change_encoding(a, bnp.encodings.DNAEncoding)

    assert np.all(a.encoding.decode(a) == b.encoding.decode(b))


def test_change_encoding_on_encoded_ragged_array():
    a = bionumpy.encoded_array_functions.as_encoded_array(["ACTG", "AC"], bnp.encodings.alphabet_encoding.ACTGEncoding)
    b = bionumpy.encoded_array_functions.change_encoding(a, bnp.encodings.DNAEncoding)

    assert_raggedarray_equal(a.encoding.decode(a.ravel()), b.encoding.decode(b.ravel()))


def test_as_encoded_on_already_encoded():
    a = bionumpy.encoded_array_functions.as_encoded_array(["ACTG"], bnp.encodings.DNAEncoding)


class CustomEncoding(OneToOneEncoding):
    def _encode(self, data):
        return data + 1

    def _decode(self, data):
        return data - 1


class CustomNumericEncoding(NumericEncoding):
    def _decode(self, data):
        return data

    def _encode(self, data):
        return data


@pytest.mark.parametrize("data", ["test", ["test1", "test2"]])
def test_public_encode_decode_string(data):
    custom_encoding = CustomEncoding()
    encoded = custom_encoding.encode(data)
    decoded = custom_encoding.decode(encoded)
    encoded2 = custom_encoding.encode(decoded)
    assert np.all(encoded == encoded2)


@pytest.mark.parametrize("data", [np.array([1, 2, 3]), np.array([[1, 10], [1]])])
def test_encode_numeric(data):
    custom_encoding = CustomNumericEncoding()
    encoded = custom_encoding.encode(data)
    decoded = custom_encoding.decode(encoded)
    encoded2 = custom_encoding.encode(decoded)
    assert np.all(encoded == encoded2)


@pytest.mark.parametrize("data",
                         ["1234",
                         ["1234", "5678"],
                         np.array([1, 2, 3, 4]),
                         RaggedArray([[1, 2, 3], [4]])])
def test_digit_encoding(data):
    encoding = DigitEncoding
    encoded = encoding.encode(data)
    decoded = encoding.decode(encoded)
    encoded2 = encoding.encode(decoded)
    assert_raggedarray_equal(encoded, encoded2)


def test1():
    seq = as_encoded_array(["ACTG", "AC"])
    print(seq)
    for s in seq:
        print("SEQ: ", s.to_string())
    assert True
