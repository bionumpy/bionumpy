import numpy as np
import pytest
from npstructures import npdataclass
from bionumpy.code_nodes import ComputationNode, LeafNode, consume, Mockup


@npdataclass
class ExampleDatatype:
    a: int
    b: int
    c: int


all_data = ExampleDatatype(np.array([0, 1, 2, 3, 4]), np.array([2, 3, 4, 1, 0]), np.array([100, 200, 1, 2, 3]))


@pytest.fixture
def whole_data():
    return all_data


@pytest.fixture
def stream():
    return Mockup((c for c in (all_data[:2], all_data[2:4], all_data[4:])))


def add(a, b, c):
    return a+b


def two_step(a, b, c):
    return a*b+b-a


def complicated(a, b, c):
    d = a+b
    e = b * c
    f = c-a > 0
    return np.where(f, np.exp(d+e+f), 1000)


@pytest.mark.parametrize("func", [add, two_step, complicated])
def test_equal(whole_data, stream, func):
    a = stream.a
    b = stream.b
    assert isinstance(a, LeafNode) and isinstance(b, LeafNode)
    true = func(whole_data.a, whole_data.b, whole_data.c)
    temp = func(a, b, stream.c)
    ours = np.concatenate(consume(temp))
    assert np.all(true == ours)
