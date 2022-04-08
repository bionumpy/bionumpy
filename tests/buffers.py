import pytest

import numpy as np

def chunk_from_text(text):
    return np.frombuffer(bytes(text, encoding="utf8"), dtype=np.uint8)

@pytest.fixture
def fastq_buffer():
    t = """\
@headerishere
CTTGTTGA
+
!!!!!!!!
@anotherheader
CGG
+
!!!
"""
    return chunk_from_text(t)


@pytest.fixture
def twoline_fasta_buffer():
    t = """\
>header
CTTGTTGA
>header2
CGG
"""
    return chunk_from_text(t)

@pytest.fixture
def bed_buffer():
    t = """\
chr1\t1\t3\t.\t.\t-
chr1\t40\t60\t.\t.\t+
chr2\t400\t600\t.\t.\t+
"""
    return chunk_from_text(t)

@pytest.fixture
def vcf_buffer():
    t = """\
chr1	88362	rs4970378	A	G
chr1	887560	rs3748595	A	C
chr2	8878	rs3828047	A	G
"""
    return chunk_from_text(t)

@pytest.fixture
def vcf_buffer2():
    t = """\
chr1	88362	rs4970378	A	G
chr1	887560	rs3748595	A	CAA
chr2	8878	rs3828047	AGG	C
"""
    return chunk_from_text(t)

@pytest.fixture
def vcf_matrix_buffer():
    return chunk_from_text("""\
chr1	883625	rs4970378	A	G\t.\t.\t.\t.\t1|1:0,4:4:6:70,6,0	1|1:0,19:19:36:358,36,0	1|1:0,3:3:6:67,6,0	1|1:0,1:1:3:34,3,0
chr1	887560	rs3748595	A	C\t.\t.\t.\t.\t0/0:7,0:7:15:0,15,163	1/1:0,30:30:81:888,81,0	1/1:0,2:2:6:68,6,0	1/1:0,1:1:3:36,3,0
chr1	887801	rs3828047	A	G\t.\t.\t.\t.\t./.	1/1:0,17:17:39:398,39,0	1/1:0,3:3:9:102,9,0	1/1:0,1:1:3:34,3,0
""")

                           
                    
