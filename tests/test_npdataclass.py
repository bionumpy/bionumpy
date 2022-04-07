import pytest
import numpy as np
import dataclasses
from bionumpy.npdataclass import NpDataClass, VarLenArray


@dataclasses.dataclass
class DummyClass(NpDataClass):
    data1: np.ndarray
    data2: np.ndarray

    def __eq__(self, other):
        return all(np.all(s == o) for s, o in zip(self.shallow_tuple(), other.shallow_tuple()))


@pytest.fixture
def objects():
    return [DummyClass(np.arange(4), 2*np.arange(4)),
            DummyClass(np.arange(3)+2, np.arange(3))]


def test_getitem(objects):
    assert objects[0][2] == DummyClass(2, 4)
    assert objects[0][:2] == DummyClass([0, 1], [0, 2])


def test_concat(objects):
    assert np.concatenate(objects) == DummyClass(np.concatenate((np.arange(4), np.arange(3)+2)),
                                                 np.concatenate((2*np.arange(4), np.arange(3))))


def test_varlenarray():
    a1 = VarLenArray(np.arange(6).reshape(3, 2))
    a2 = VarLenArray(np.arange(3).reshape(3, 1))
    v = np.concatenate((a1, a2))
    assert np.all(v.array == np.concatenate((a1.array, [[0, 0], [0, 1], [0, 2]])))
