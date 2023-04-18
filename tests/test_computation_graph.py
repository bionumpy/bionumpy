import numpy as np
from numpy.testing import assert_equal
import pytest
from npstructures import npdataclass
from bionumpy.computation_graph import StreamNode, compute
# from bionumpy.code_nodes import LeafNode, consume, NpDataclassStream


@npdataclass
class ExampleDatatype:
    a: int
    b: int
    c: int


@pytest.fixture
def data():
    return (np.array([0, 1, 2, 3, 4]), np.array([2, 3, 4, 1, 0]), np.array([100, 200, 1, 2, 3]))


@pytest.fixture
def stream_data(data):
    return tuple(StreamNode(iter([col[:2], col[2:4], col[4:]])) for col in data)


@pytest.fixture
def mixed_data(data, stream_data):
    return (stream_data[0], data[1], stream_data[2])


def add(a, b, c):
    return a+b


def two_step(a, b, c):
    return a*b+b-a


def complicated(a, b, c):
    d = a+b
    print(str(d))
    e = b * c
    print(str(e))
    f = c-a > 0
    print(str(f))
    g = np.where(f, np.exp(d+e+f), 1000)
    print(g)
    return(g)


@pytest.mark.parametrize("func", [add, two_step, complicated])
def test_equal(data, stream_data, func):
    true = func(*data)
    ours = func(*stream_data).compute()
    assert_equal(true, ours)


@pytest.mark.parametrize("func", [add, two_step, complicated])
@pytest.mark.skip('should fail')
def test_mixed_equal(data, mixed_data, func):
    true = func(*data)
    ours = func(*mixed_data).compute()
    assert_equal(true, ours)


def test_print(stream_data):
    a, b, c = stream_data
    print(a, b)
    print(complicated(a, b, c))


def test_reduce(stream_data, data):
    stream_data = stream_data[0]
    data = data[0]
    true = np.sum(data)
    s = np.sum(stream_data)
    solution = s.compute()
    assert true == solution


def test_histogram(stream_data, data):
    stream_data = stream_data[0]
    data = data[0]
    true = np.histogram(data, bins=3, range=(0, 6))
    s = np.histogram(stream_data, bins=3, range=(0, 6))
    solution = s.compute()
    assert_equal(true[0], solution[0])

@pytest.mark.xfail
def test_double_compute(stream_data, data):
    stream_data = stream_data[0]
    a = stream_data+1
    b = stream_data+2
    c = b+2
    a_true = data[0]+1
    b_true = data[0]+2
    c_true = b_true+2
    a, b, c = compute([a, b, c])
    assert_equal(a, a_true)
    assert_equal(b, b_true)
    assert_equal(c, c_true)


def test_double_reduce(stream_data, data):
    stream_data = stream_data[0]
    data = data[0]
    true_hist = np.histogram(data, bins=3, range=(0, 6))
    true_sum = np.sum(data)
    h = np.histogram(stream_data, bins=3, range=(0, 6))
    s = np.sum(stream_data)
    h, s = compute((h, s))
    assert_equal(h[0], true_hist[0])
    assert_equal(s, true_sum)
