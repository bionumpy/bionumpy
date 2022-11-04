import pytest
from bionumpy.io.delimited_buffers import BedBuffer
from bionumpy.datatypes import Interval
intervals = Interval(["chr1", "chr1"], [2, 10], [100, 20])


def test_from_data():
    buf = BedBuffer.from_data(intervals)
    assert buf.to_string() == """\
chr1\t2\t100
chr1\t10\t20
"""
