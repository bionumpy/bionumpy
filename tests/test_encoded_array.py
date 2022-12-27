import pytest
from bionumpy.encoded_array import EncodedArray, EncodedRaggedArray
from bionumpy import as_encoded_array
from bionumpy.encodings import BaseEncoding, QualityEncoding
from bionumpy.util.testing import assert_encoded_raggedarray_equal, assert_encoded_array_equal
import bionumpy as bnp
import numpy as np


@pytest.fixture
def sequence():
    return EncodedArray(np.array([[1, 2], [2, 3]], dtype=np.uint8)+64, bnp.encodings.BaseEncoding)

@pytest.fixture
def simple_sequence():
    return EncodedArray(np.array([1, 2, 2, 3], dtype=np.uint8) + 64, bnp.encodings.BaseEncoding)


@pytest.fixture
def ragged_string_list():
    return ['a', 'ac', 'acg']


def _test_sequence_repr(sequence):
    assert repr(sequence).replace("\n", "") == "EncodedArray([[65, 66] [66 67]])"


def test_sequence_str(sequence):
    assert str(sequence) == """['AB' 'BC']"""


def test_concat_sequenc(sequence):
    cat = np.concatenate([sequence, sequence])
    assert isinstance(cat, EncodedArray)


def test_iter(sequence):
    for elem in sequence.ravel():
        assert isinstance(elem, EncodedArray)


def test_getitem(sequence):
    sequence = sequence.ravel()
    for i in range(len(sequence)):
        assert isinstance(sequence[i], EncodedArray)


def test_works_with_ragged():
    encoded_array = EncodedArray([67, 68, 69], bnp.encodings.BaseEncoding)
    encoded_array == "D"
    ragged = EncodedRaggedArray(encoded_array, [2, 1])


def test_empty_encoded():
    result = as_encoded_array(['', ''])
    true = EncodedRaggedArray(EncodedArray(np.empty((0,), dtype=np.uint8), BaseEncoding),
                              [0, 0])
                                           
    assert_encoded_raggedarray_equal(result, true)


def test_empty_dna_encoded():
    result = as_encoded_array(['', ''], bnp.DNAEncoding)
    true = EncodedRaggedArray(EncodedArray(np.empty((0,), dtype=np.uint8), bnp.DNAEncoding),
                              [0, 0])
                                           
    assert_encoded_raggedarray_equal(result, true)

def test_array_function(simple_sequence):
    result = np.append(simple_sequence, simple_sequence)
    concat_sequence = EncodedArray(np.array([1, 2, 2, 3, 1, 2, 2, 3], dtype=np.uint8) + 64, bnp.encodings.BaseEncoding)
    assert_encoded_array_equal(result,  concat_sequence)

def test_insert_function(simple_sequence):
    result = np.insert(simple_sequence, 0, simple_sequence)
    concat_sequence = EncodedArray(np.array([1, 2, 2, 3, 1, 2, 2, 3], dtype=np.uint8) + 64, bnp.encodings.BaseEncoding)
    assert_encoded_array_equal(result, concat_sequence)


def test_encoded_array_list_to_raggedarray(simple_sequence):
    array_list = as_encoded_array([simple_sequence, simple_sequence])
    str_list = as_encoded_array([simple_sequence.to_string(), simple_sequence.to_string()])
    assert_encoded_raggedarray_equal(array_list, str_list)
                                

@pytest.mark.parametrize('dtype', [object, None])
def test_object_array(ragged_string_list, dtype):
    object_array = np.array(ragged_string_list, dtype=dtype)
    encoded_array = as_encoded_array(object_array)
    assert_encoded_raggedarray_equal(encoded_array, ragged_string_list)
    

if __name__ == "__main__":
    test_works_with_ragged()
