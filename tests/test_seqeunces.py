import pytest
from bionumpy.sequences import Sequence
import numpy as np

@pytest.fixture
def sequence():
    return (np.array([[1, 2], [2, 3]], dtype=np.uint8)+64).view(Sequence)

def test_sequence_repr(sequence):
    
    assert repr(sequence).replace("\n", "") == "Sequence([[65 66] [66 67]], ASCIIEncoding)"


def test_sequence_str(sequence):
    assert str(sequence) == """['AB' 'BC']"""
