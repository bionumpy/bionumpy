import numpy as np
from npstructures.testing import assert_raggedarray_equal

import bionumpy as bnp
from tests.io_tests.sub_delimiters import find_sub_delimiters


def test_find_sub_delimiters():
    text = bnp.as_encoded_array('headerhere\ta:b:c\notherheader\td:e\n')
    starts = np.flatnonzero(text=='\t')+1
    ends = np.flatnonzero(text=='\n')
    result = find_sub_delimiters(text, starts, ends, delimiter=':')
    assert_raggedarray_equal(result, [starts[0]+np.array([1, 3]),
                                       starts[1]+np.array([1])])

