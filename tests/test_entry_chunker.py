from bionumpy.streams.chunk_entries import chunk_entries
from bionumpy.streams import BnpStream
import numpy as np
import pytest


@pytest.fixture
def array_chunks():
    return BnpStream([np.arange(4), np.arange(3), np.arange(2), np.arange(1)])


def test_chunk_entries(array_chunks):
    chunked = list(chunk_entries(array_chunks, 3))
    assert all(len(chunk) == 3 for chunk in chunked[:-1])
    assert len(chunked[-1]) <= 3
    np.testing.assert_array_equal(
        np.concatenate(chunked),
        np.concatenate([np.arange(4), np.arange(3), np.arange(2), np.arange(1)]))
