import dataclasses
from numpy.testing import assert_array_equal, assert_array_almost_equal
from npstructures.testing import assert_raggedarray_equal
from npstructures.npdataclasses import shallow_tuple
from ..encoded_array import EncodedArray, EncodedRaggedArray, as_encoded_array
from npstructures import RaggedArray
import numpy as np

def assert_encoded_array_equal(array1, array2):
    array1, array2 = (as_encoded_array(a) for a in (array1, array2))
    assert array1.encoding == array2.encoding, (array1.encoding, array2.encoding)
    assert_array_equal(array1.raw(), array2.raw())

def assert_encoded_raggedarray_equal(array1, array2):
    array1, array2 = (as_encoded_array(a) for a in (array1, array2))
    assert array1.encoding == array2.encoding, (array1.encoding, array2.encoding)
    assert_raggedarray_equal(array1.raw(), array2.raw())

def assert_float_close_enough(a, b):
    if np.allclose(a, b):
        return
    fa, ma = np.frexp(a)
    fb, mb = np.frexp(b)
    fa = np.where(ma>mb, fa*2**(np.maximum(ma-mb, 0)), fa)
    fb = np.where(mb>ma, fb*2**(np.maximum(mb-ma, 0)), fb)
    assert_array_almost_equal(fa, fb)

def assert_bnpdataclass_equal(a, b):
    assert dataclasses.fields(a) == dataclasses.fields(b)
    for s, o, field in zip(shallow_tuple(a), shallow_tuple(b), dataclasses.fields(a)):
        assert issubclass(type(s), type(o)) or issubclass(type(o), type(s)), (type(s), type(o))
        if isinstance(s, EncodedRaggedArray):
            assert_encoded_raggedarray_equal(s, o)
        elif isinstance(s, EncodedArray):
            assert_encoded_array_equal(s, o)
        elif isinstance(s, RaggedArray):
            assert_raggedarray_equal(s, o)
        elif isinstance(s, np.ndarray):
            if field.type == float:
                assert_float_close_enough(s, o)
            else:
                assert_array_equal(s, o)
        else:
            assert hasattr(s, "shape") and hasattr(o, "shape"), (s, o)
            if not s.shape == o.shape:
                assert False, (s.shape, o.shape, str(a), str(b))
            if not np.all(np.equal(s, o)):
                assert False, (a, b, field.name)
