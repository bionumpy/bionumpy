import gzip
import numpy as np
from .buffers import data
import pytest
from bionumpy.bam import BamBuffer, alignment_to_interval


@pytest.mark.skip("bam")
def test_read():
    filename = "/home/knut/Sources/bionumpy/example_data/test.bam"
    f = gzip.open(filename, "rb")
    header_data = BamBuffer.read_header(f)
    chunk = np.frombuffer(f.read(3000), dtype=np.uint8)
    buf = BamBuffer.from_raw_buffer(chunk, header_data=header_data)
    d = buf.get_data()
    print(d)
    print(alignment_to_interval(d))
    for a, b in zip(data["sam"], d):
        assert b == a, (b, a)
        # assert a.position == b.position
        # assert a.name == b.name
    assert False
