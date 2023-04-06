import pytest
from bionumpy.sequence.debruin import DeBruijnGraph
from bionumpy.util.testing import assert_encoded_array_equal

@pytest.fixture
def sequences():
    return ['acg', 'gtc']

@pytest.fixture
def sequences():
    return ['acg', 'gtc']




def test_create_debruijn(sequences):
    graph = DeBruijnGraph.from_sequences(sequences, 2)
    assert isinstance(graph, DeBruijnGraph)


def test_forward(sequences):
    graph = DeBruijnGraph.from_sequences(sequences, 2)
    assert graph.forward('ac') == ['CG']


def test_backward(sequences):
    graph = DeBruijnGraph.from_sequences(sequences, 2)
    assert graph.backward('tc') == ['GT']


@pytest.mark.xfail
def test_compatibility_classes():
    pass


@pytest.mark.xfail
def test_():
    pass


def test_fw(sequences):
    pass

