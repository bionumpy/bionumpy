import pytest
import numpy as np
from bionumpy.groupby import get_changes
from npstructures import RaggedArray


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
        [[2, 3, 4]]]


def test_get_changes(ragged_array):
    changes = get_changes(ragged_array)
    np.testing.assert_equal(changes, [2, 4, 6])
