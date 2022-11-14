from bionumpy.streams import NpDataclassStream, quantile, mean
import numpy as np
import dataclasses
import pytest


@dataclasses.dataclass
class DummyClass:
    pos: int


@pytest.fixture
def test_stream():
    return NpDataclassStream(DummyClass(word) for word in ["hei", "pa", "deg"])


def test_str(test_stream):
    s = str(test_stream)
    assert s.endswith("DummyClass(pos='hei')")
    assert "DummyClass" in s


def test_quantile():
    array = [1, 5, 10]
    assert quantile(array, 0.5) == 5


def test_quantile1():
    array = [1, 10]
    assert quantile(array, 0.5) == 1


def test_quantile2():
    array = [1, 10, 12, 20]
    assert quantile(array, 0.5) == 10


def test_mean():
    array = np.array([1, 5, 10])
    assert mean(array) == 16/3


def test_mean1():
    array = np.array([1, 10])
    assert mean(array) == 5.5


def test_mean2():
    array = np.array([1, 10, 12, 20])
    assert mean(array) == 43/4


@pytest.mark.skip("Unimplemented")
def test_colmean():
    array = np.array([[1, 10], [12, 20]])
    np.testing.assert_array_equal(mean(array, axis=0), [5.5, 16])
