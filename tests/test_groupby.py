import pytest
import numpy as np
from bionumpy.streams.groupby_func import get_changes, groupby, join_groupbys
from bionumpy.encoded_array import EncodedArray, EncodedRaggedArray
from npstructures import RaggedArray
from npstructures.testing import assert_raggedarray_equal
import bionumpy as bnp

@pytest.fixture
def ragged_array():
    return RaggedArray([
        [1, 2, 3],
        [1, 2, 3],
        [1, 2],
        [1, 2],
        [2, 3, 4],
        [2, 3, 4],
        [2, 3, 3]])


@pytest.fixture
def grouped():
    return [
        [[1, 2, 3],
         [1, 2, 3]],
        [[1, 2],
         [1, 2]],
        [[2, 3, 4],
         [2, 3, 4]],
        [[2, 3, 3]]]


@pytest.mark.parametrize("do_encode", [True, False])
def test_get_changes(ragged_array, do_encode):
    if do_encode:
        ragged_array = EncodedRaggedArray(EncodedArray(ragged_array.ravel(), bnp.encodings.BaseEncoding), ragged_array.shape)
    changes = get_changes(ragged_array)
    np.testing.assert_equal(changes, [2, 4, 6])

@pytest.mark.parametrize("do_encode", [True, False])
def test_groupby(ragged_array, grouped, do_encode):
    if do_encode:
        ragged_array = EncodedRaggedArray(EncodedArray(ragged_array.ravel(), bnp.encodings.BaseEncoding), ragged_array.shape)
    groups = list(groupby(ragged_array))
    assert len(groups) == len(grouped)
    for (k, g), true in zip(groups, grouped):
        true = RaggedArray(true)
        if do_encode:
            true = EncodedRaggedArray(EncodedArray(true.ravel(), bnp.encodings.BaseEncoding), true.shape)
        assert_raggedarray_equal(g, true)

@pytest.mark.parametrize("do_encode", [True, False])
def test_join_groupbys(ragged_array, grouped, do_encode):
    if do_encode:
        ragged_array = EncodedRaggedArray(EncodedArray(ragged_array.ravel(), bnp.encodings.BaseEncoding), ragged_array.shape)
    gs = [groupby(ragged_array), groupby(ragged_array)]
    groups = list(join_groupbys(gs))
    print(groups)
    for (k, g), true in zip(groups, grouped*2):
        true = RaggedArray(true)
        if do_encode:
            true = EncodedRaggedArray(EncodedArray(true.ravel(), bnp.encodings.BaseEncoding), true.shape)

        assert_raggedarray_equal(g, true)
