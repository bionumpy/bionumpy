import pytest
from bionumpy.sequence.debruin import DeBruijnGraph, ColoredDeBruijnGraph
from bionumpy.util.testing import assert_encoded_array_equal


@pytest.fixture
def sequences():
    return ['acg', 'cgtc']


def test_create_debruijn(sequences):
    graph = DeBruijnGraph.from_sequences(sequences, 2)
    assert isinstance(graph, DeBruijnGraph)


def test_forward(sequences):
    graph = DeBruijnGraph.from_sequences(sequences, 2)
    assert graph.forward('ac') == ['CG']


def test_backward(sequences):
    graph = DeBruijnGraph.from_sequences(sequences, 2)
    assert graph.backward('tc') == ['GT']


def test_create_colored_debruijn_graph(sequences):
    graph = ColoredDeBruijnGraph.from_sequences(sequences, 2)
    assert isinstance(graph, ColoredDeBruijnGraph)


def test_colored_debruijn_graph(sequences):
    graph = ColoredDeBruijnGraph.from_sequences(sequences, 2)
    assert graph['ac'] == [0]
    assert graph['tc'] == [1]
    assert graph['cg'] == [0, 1]



@pytest.mark.xfail
def test_compatibility_classes():
    pass


@pytest.mark.xfail
def test_():
    pass


def test_fw(sequences):
    pass

