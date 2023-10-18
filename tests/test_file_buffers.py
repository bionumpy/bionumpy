import numpy as np
import pytest

from bionumpy import as_encoded_array
from bionumpy.io.file_buffers import TextThroughputExtractor
from bionumpy.util.testing import assert_encoded_array_equal


@pytest.fixture
def data():
    return as_encoded_array('hei,pa,deg')

@pytest.fixture
def buffer_extractor(data):
    return TextThroughputExtractor(data,
                                   field_starts = np.array([0, 4, 7])[:, None],
                                   field_ends = np.array([3, 6, 10])[:, None],
                                   entry_starts = np.array([0, 4, 7]),
                                   entry_ends = np.array([4, 7, 10]))


def test_data(buffer_extractor):
    subset = buffer_extractor[[0, 2]]
    assert_encoded_array_equal(subset.data, 'hei,deg')



