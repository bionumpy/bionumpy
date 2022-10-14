import pytest
from bionumpy.encoded_array import EncodedArray, EncodedRaggedArray
import bionumpy as bnp




def test_works_with_ragged():
    encoded_array = EncodedArray([67, 68, 69], bnp.encodings.BaseEncoding)
    encoded_array == "D"
    ragged = EncodedRaggedArray(encoded_array, [2, 1])
    assert False

if __name__ == "__main__":
    test_works_with_ragged()
