import pytest
import gzip
import numpy as np
from .buffers import data


from bionumpy.bam import BamBuffer, BAMReader


def test_read():
    filename = "/home/knut/Sources/bionumpy/example_data/alignments.bam"
    f = gzip.open(filename, "rb")
    BAMReader(f)
    chunk = np.frombuffer(f.read(300), dtype=np.uint8)
    buf = BamBuffer.from_raw_buffer(chunk)
    d = buf.get_data()
    for a, b in zip(data["sam"], d):
        assert a.position == b.position
