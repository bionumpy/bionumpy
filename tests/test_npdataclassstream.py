from bionumpy.npdataclassstream import NpDataclassStream
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
