import pytest
import numpy as np
from bionumpy.bam import BamBuffer, BAMReader


def test_read():
    filename = "/home/knut/Downloads/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam"
    f = open(filename, "rb")
    BAMReader(f)
    chunk = np.frombuffer(f.read(300), dtype=np.uint8)
    buf = BamBuffer.from_raw_buffer(chunk)
    buf.get_data()
    assert False
