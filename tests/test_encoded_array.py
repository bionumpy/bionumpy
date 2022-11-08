import pytest
from bionumpy.encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from bionumpy.testing import assert_encoded_raggedarray_equal
import bionumpy as bnp
import numpy as np


@pytest.fixture
def sequence():
    return EncodedArray(np.array([[1, 2], [2, 3]], dtype=np.uint8)+64, bnp.encodings.BaseEncoding)


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
    true = EncodedRaggedArray(EncodedArray(np.empty((0,), dtype=np.uint8)),
                              [0, 0])
                                           
    assert_encoded_raggedarray_equal(result, true)


def test_empty_dna_encoded():
    result = as_encoded_array(['', ''], bnp.DNAEncoding)
    true = EncodedRaggedArray(EncodedArray(np.empty((0,), dtype=np.uint8), bnp.DNAEncoding),
                              [0, 0])
                                           
    assert_encoded_raggedarray_equal(result, true)
    


if __name__ == "__main__":
    test_works_with_ragged()
