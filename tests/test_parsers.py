import numpy as np
from bionumpy.parser import FastQBuffer

def chunk_from_text(text):
    return np.frombuffer(bytes(text, encoding="utf8"), dtype=np.uint8)

def test_fastq_buffer():
    t = """@headerishere
CTTGTTGA
+
!!!!!!!!
@anotherheader
CGG
+
!!!
"""
    chunk = chunk_from_text(t)
    buf = FastQBuffer.from_raw_buffer(chunk)
    seqs = buf.get_sequences()
    assert [buf._encoding.to_string(seq) for seq in seqs] == ["CTTGTTGA", "CGG"]
