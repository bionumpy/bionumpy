import bionumpy as bnp
import pytest
from bionumpy.io.bam import BamIntervalBuffer
from bionumpy.alignments import alignment_to_interval


#@pytest.mark.skip
def test_read():
    filename = "example_data/test.bam"
    f = bnp.open(filename)
    d = f.read()
    print(d)
    print(d.flag.dtype)


def test_read_intervals():
    filename = "example_data/test.bam"
    f = bnp.open(filename, buffer_type=BamIntervalBuffer)
    d = f.read()
    print(d)
