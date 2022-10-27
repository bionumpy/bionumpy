from numpy.testing import assert_array_equal
from npstructures.testing import assert_raggedarray_equal

def assert_encoded_array_equal(array1, array2):
    assert array1.encoding == array2.encoding, (array1.encoding, array2.encoding)
    assert_array_equal(array1.raw(), array2.raw())

def assert_encoded_raggedarray_equal(array1, array2):
    assert array1.encoding == array2.encoding, (array1.encoding, array2.encoding)
    assert_raggedarray_equal(array1.raw(), array2.raw())
