from bionumpy.npdataclassstream import NpDataclassStream, quantile
import dataclasses
import pytest

@dataclasses.dataclass
class TestClass:
    pos: int

@pytest.fixture
def test_stream():
    return NpDataclassStream(iter(["hei", "pa", "deg"]), TestClass)


def test_str(test_stream):
    s = str(test_stream)
    assert s.endswith("hei")
    assert "TestClass" in s


def test_quantile():
    array = [1, 5, 10]
    assert quantile(array, 0.5) == 5


def test_quantile1():
    array = [1, 10]
    assert quantile(array, 0.5) == 1


def test_quantile2():
    array = [1, 10, 12, 20]
    assert quantile(array, 0.5) == 10
