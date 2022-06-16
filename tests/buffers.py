import pytest
from bionumpy.file_buffers import FastQBuffer, TwoLineFastaBuffer
from bionumpy.datatypes import SequenceEntry, SequenceEntryWithQuality, Interval, SNP
from bionumpy.delimited_buffers import BedBuffer, VCFBuffer, GfaSequenceBuffer

import numpy as np


def chunk_from_text(text):
    return np.frombuffer(bytes(text, encoding="utf8"), dtype=np.uint8)


buffer_texts = {
    "fastq": """\
@headerishere
CTTGTTGA
+
!!!!!!!!
@anotherheader
CGG
+
~~~
"""
    , "fasta": """\
>header
CTTGTTGA
>header2
CGG
"""
    , "bed": """\
chr1\t1\t3\t.\t.\t-
chr1\t40\t60\t.\t.\t+
chr2\t400\t600\t.\t.\t+
"""
    , "vcf": """\
chr1	88362	rs4970378	A	G
chr1	887560	rs3748595	A	C
chr2	8878	rs3828047	A	G
"""
    , "vcf2": """\
chr1	88362	rs4970378	A	G
chr1	887560	rs3748595	A	CAA
chr2	8878	rs3828047	AGG	C
"""
    , "vcf_matrix": """\
chr1	883625	rs4970378	A	G\t.\t.\t.\t.\t1|1:0,4:4:6:70,6,0	1|1:0,19:19:36:358,36,0	1|1:0,3:3:6:67,6,0	1|1:0,1:1:3:34,3,0
chr1	887560	rs3748595	A	C\t.\t.\t.\t.\t0/0:7,0:7:15:0,15,163	1/1:0,30:30:81:888,81,0	1/1:0,2:2:6:68,6,0	1/1:0,1:1:3:36,3,0
chr1	887801	rs3828047	A	G\t.\t.\t.\t.\t./.	1/1:0,17:17:39:398,39,0	1/1:0,3:3:9:102,9,0	1/1:0,1:1:3:34,3,0
"""
    , "gfa_sequence": """\
S\tid1\tAACCTTGG\t.\t.
S\tid4\tACTG\t*\t*
""", "gff_entry": """\
CHROMOSOME_I	Allele	substitution	10017380	10017380	.	+.	aachange=A to T;consequence=Missense;interpolated_map_position=4.49151;public_name=q504;substitution=G/A;variation=WBVar00241143
CHROMOSOME_I	Allele	substitution	10573196	10573196	.	+.	aachange=A to T;consequence=Missense;interpolated_map_position=5.05579;public_name=vc56;substitution=G/A;variation=WBVar00275020
CHROMOSOME_I	Allele	substitution	13447169	13447169	.	+.	aachange=A to T;consequence=Missense;interpolated_map_position=20.7564;public_name=ar197;strain=GS3412;substitution=G/A;variation=WBVar00000227
CHROMOSOME_I	Allele	substitution	13721389	13721389	.	+.	aachange=A to V;consequence=Missense;interpolated_map_position=21.6433;public_name=mg334;substitution=G/A;variation=WBVar00088937
"""
}

buffers = {key: chunk_from_text(val) for key, val in buffer_texts.items()}

data = {
        "bed": [
            Interval("chr1", 1, 3),
            Interval("chr1", 40, 60),
            Interval("chr2",  400, 600)],
        "vcf2": [
            SNP("chr1",	88361,	"A",	"G"),
            SNP("chr1",	887559,	"A",	"CAA"),
            SNP("chr2",	8877,	"AGG",	"C")],
        "vcf": [
            SNP("chr1",	88361,	"A",	"G"),
            SNP("chr1",	887559,	"A",	"C"),
            SNP("chr2",	8877,	"A",	"G")],
        "fastq": [
            SequenceEntryWithQuality("headerishere", "CTTGTTGA", [0 for _ in "CTTGTTGA"]),
            SequenceEntryWithQuality("anotherheader", "CGG", [ord("~")-ord("!") for _ in "CGG"])],
        "fasta": [
            SequenceEntry("header", "CTTGTTGA"),
            SequenceEntry("header2", "CGG")],
        "gfa_sequence": [
            SequenceEntry("id1", "AACCTTGG"),
            SequenceEntry("id4", "ACTG")
        ]
}


buffer_type = {"bed": BedBuffer,
               "vcf2": VCFBuffer,
               "vcf": VCFBuffer,
               "fastq": FastQBuffer,
               "fasta": TwoLineFastaBuffer,
               "gfa_sequence": GfaSequenceBuffer}

combos = {key: (buffers[key], data[key], buffer_type[key]) for key in buffer_type}



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

@pytest.fixture
def gfa_sequence_buffer():
    t = """\
S\tid1\tAACCTTGG\t.\t.
S\tid4\tACTG\t*\t*
"""
    return chunk_from_text(t)
