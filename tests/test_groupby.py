import pytest
import numpy as np
from bionumpy.groupby import get_changes, groupby, join_groupbys
from npstructures import RaggedArray
from npstructures.testing import assert_raggedarray_equal


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


def test_get_changes(ragged_array):
    changes = get_changes(ragged_array)
    np.testing.assert_equal(changes, [2, 4, 6])


def test_groupby(ragged_array, grouped):
    groups = list(groupby(ragged_array))
    assert len(groups) == len(grouped)
    for (k, g), true in zip(groups, grouped):
        assert_raggedarray_equal(g, true)


def test_join_groupbys(ragged_array, grouped):
    gs = [groupby(ragged_array), groupby(ragged_array)]
    groups = list(join_groupbys(gs))
    print(groups)
    for (k, g), true in zip(groups, grouped*2):
        assert_raggedarray_equal(g, true)
