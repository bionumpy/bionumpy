import bionumpy as bnp
import pytest
from bionumpy.io.bam import BamIntervalBuffer
from bionumpy.alignments import alignment_to_interval


#@pytest.mark.skip
def test_read_acceptance():
    filename = "example_data/test.bam"
    f = bnp.open(filename)
    d = f.read()
    print(d)
    print(d.flag.dtype)


def test_read_intervals_acceptance():
    filename = "example_data/test.bam"
    f = bnp.open(filename, buffer_type=BamIntervalBuffer)
    d = f.read()
    print(d)

@pytest.mark.xfail
def test_read_bam():
    filename = 'example_data/small_alignments.bam'
    entries = bnp.open(filename).read()
    assert entries[0].position == 3837783
