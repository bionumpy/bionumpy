import bionumpy as bnp
import gzip
import numpy as np
from .buffers import data
import pytest
from bionumpy.bam import BamBuffer, alignment_to_interval


@pytest.mark.skip
def test_read():
    filename = "/home/knut/Sources/bionumpy/example_data/test.bam"
    f = bnp.open(filename)
    d = f.read()
    print(d)
    # d = np.concatenate(list(f.read_chunks(1000)))
    print(alignment_to_interval(d))
    # assert False
