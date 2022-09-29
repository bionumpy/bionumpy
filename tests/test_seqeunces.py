import pytest
from bionumpy.sequences import Sequence
import numpy as np

@pytest.fixture
def sequence():
    return (np.array([[1, 2], [2, 3]], dtype=np.uint8)+64).view(Sequence)

def _test_sequence_repr(sequence):
    assert repr(sequence).replace("\n", "") == "Sequence([[65, 66] [66 67]])"


def test_sequence_str(sequence):
    assert str(sequence) == """['AB' 'BC']"""

def test_concat_sequenc(sequence):
    cat = np.concatenate([sequence, sequence])
    assert isinstance(cat, Sequence)


def test_iter(sequence):
    for elem in sequence.ravel():
        assert isinstance(elem, Sequence), (elem, type(elem), elem.view(Sequence))


def test_getitem(sequence):
    sequence = sequence.ravel()
    for i in range(len(sequence)):
        assert isinstance(sequence[i], Sequence)
