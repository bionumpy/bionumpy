from bionumpy.util.ragged_slice import ragged_slice
import bionumpy as bnp
import numpy as np

def test_ragged_slice():
    data = bnp.as_encoded_array(["ACTG", "AAC", "A"])
    sliced = ragged_slice(data, starts=np.array([1, 0, 0]), ends=np.array([3, 2, 1]))

    assert np.all(sliced == bnp.as_encoded_array(["CT", "AC", "A"]))
