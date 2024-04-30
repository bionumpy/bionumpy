import bionumpy as bnp
from tempfile import NamedTemporaryFile
import numpy as np


def test_read_write_pairs_file(data_path):
    file = data_path / "small.pairs"
    data = bnp.open(file).read()

    assert data[0].pos1 == 61

    file = NamedTemporaryFile(suffix=".pairs")
    with bnp.open(file.name, mode='w') as f:
        f.write(data)

    new_data = bnp.open(file.name).read()
    assert np.all(data == new_data)