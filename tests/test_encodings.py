import numpy as np
import pytest
import bionumpy as bnp
from npstructures.testing import assert_raggedarray_equal
from npstructures import RaggedArray
import bionumpy.encoded_array
import bionumpy.encoded_array_functions
import bionumpy as bnp
from bionumpy import FastQBuffer
from bionumpy.datatypes import SequenceEntryWithQuality
from bionumpy.encodings import DigitEncoding, QualityEncoding, CigarEncoding
from bionumpy.encodings.base_encoding import NumericEncoding, OneToOneEncoding, BaseEncoding
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


@pytest.mark.parametrize("data", ["!!@-^", ["!!@@", "!+"]])
def test_base_quality_encoding(data):
    encoding = QualityEncoding
    encoded = encoding.encode(data)
    print(type(encoded))
    print(encoded)
    decoded = encoding.decode(encoded)
    print(type(decoded))
    print(decoded)
    encoded2 = encoding.encode(decoded)
    assert np.all(encoded2 == encoded)


def test1():
    seq = as_encoded_array(["ACTG", "AC"])
    print(seq)
    for s in seq:
        print("SEQ: ", s.to_string())
    assert True


def test_encoding_sequence_entry():
    s = SequenceEntryWithQuality(
        name=as_encoded_array(['headerishere'], BaseEncoding),
        sequence=as_encoded_array(['CTTGTTGA'], BaseEncoding),
        quality=np.array([[223, 223, 223, 223, 223, 223, 223, 223]], dtype=np.uint8))

    data = FastQBuffer.dataclass.stack_with_ragged(s)
    print("DATA")
    #print(data)
    #buf = FastQBuffer.from_data(data)
    #print(buf.raw())


@pytest.mark.parametrize("data", ["!!@-^", ["!!@@", "!+"]])
def test_cigar_encoding(data):

    encoding = CigarEncoding
    encoded = encoding.encode(data)
    encoded_old = as_encoded_array(data, encoding)
    print(encoded)
    print(encoded_old)
    assert_raggedarray_equal(encoded, encoded_old)
    decoded = encoding.decode(encoded)
    print(decoded)
    encoded2 = encoding.encode(decoded)
    encoded2_old = as_encoded_array(decoded, encoding)
    #assert assert_raggedarray_equal(encoded2, encoded2_old)
    print(encoded2)
    #assert np.all(encoded2 == encoded)



class TestNumericEncoding(NumericEncoding):
    def _encode(self, data):
        return data - 49

    def _decode(self, data):
        return data + 49


def test_custom_numeric_encoding():
    encoding = QualityEncoding
    string = "!#!#"
    array = np.array([33, 35, 33, 35])
    encoded = encoding.encode(string)
    encoded_old = as_encoded_array(string, encoding)
    print("ENcoded new")
    print(repr(encoded))
    print("ENcoded old")
    print(repr(encoded_old))
    print("Encoded old array")
    print(as_encoded_array(array, encoding))
    decoded = encoding.decode(encoded)
    print(repr(decoded))
    #assert False