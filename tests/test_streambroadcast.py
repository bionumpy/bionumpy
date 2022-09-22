import numpy as np
import pytest
from npstructures import npdataclass
from bionumpy.code_nodes import ComputationNode, LeafNode, consume, Mockup


@npdataclass
class ExampleDatatype:
    a: int
    b: int


all_data = ExampleDatatype(np.array([0, 1, 2, 3, 4]), np.array([2, 3, 4, 1, 0]))


@pytest.fixture
def whole_data():
    return all_data


@pytest.fixture
def stream():
    return Mockup((c for c in (all_data[:2], all_data[2:4], all_data[4:])))


def add(a, b):
    return a+b


def test_equal(whole_data, stream):
    a = stream.a
    b = stream.b
    assert isinstance(a, LeafNode) and isinstance(b, LeafNode)
    true = add(whole_data.a, whole_data.b)
    temp = add(a, b)
    ours = np.concatenate(consume(temp))
    assert np.all(true == ours)
