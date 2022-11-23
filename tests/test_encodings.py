import numpy as np
import bionumpy as bnp
from npstructures.testing import assert_raggedarray_equal


def test_change_encoding_on_encoded_array():
    a = bnp.as_encoded_array("ACTG", bnp.encodings.alphabet_encoding.ACTGEncoding)
    b = bnp.encoded_array.change_encoding(a, bnp.encodings.DNAEncoding)

    assert np.all(a.encoding.decode(a) == b.encoding.decode(b))


def test_change_encoding_on_encoded_ragged_array():
    a = bnp.as_encoded_array(["ACTG", "AC"], bnp.encodings.alphabet_encoding.ACTGEncoding)
    b = bnp.encoded_array.change_encoding(a, bnp.encodings.DNAEncoding)

    assert_raggedarray_equal(a.encoding.decode(a.ravel()), b.encoding.decode(b.ravel()))


def test_as_encoded_on_already_encoded():
    a = bnp.as_encoded_array(["ACTG"], bnp.encodings.DNAEncoding)

