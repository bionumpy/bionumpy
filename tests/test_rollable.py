import pytest
import numpy as np
from bionumpy.sequence.rollable import RollableFunction
from bionumpy.encodings import BaseEncoding
from npstructures import RaggedArray


@pytest.fixture
def sequences():
    return ["acgt", "ggt", "gg"]


class Func(RollableFunction):
    _encoding = BaseEncoding

    @property
    def window_size(self):
        return 2

    def __call__(self, sequence):
        return np.all(sequence == "gt", axis=-1)


def test_rollable_same(sequences):
    result = Func().rolling_window(sequences, mode="same")
    ra = RaggedArray([[False, False, True, False],
                      [False, True, False],
                      [False, False]])
    assert np.all(result == ra)


def test_rollable_valid(sequences):
    ra = RaggedArray([[False, False, True],
                      [False, True],
                      [False]])
    result = Func().rolling_window(sequences, mode="valid")
    assert np.all(result == ra)
