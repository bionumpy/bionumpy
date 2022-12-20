import pytest
from bionumpy.io.file_buffers import FastQBuffer, TwoLineFastaBuffer
from bionumpy.datatypes import SequenceEntry, SequenceEntryWithQuality, Interval, SNP, SAMEntry, VCFEntry, Bed12, Bed6
from bionumpy.io.delimited_buffers import BedBuffer, VCFBuffer, GfaSequenceBuffer, Bed12Buffer, Bed6Buffer
from bionumpy.io.multiline_buffer import MultiLineFastaBuffer
from bionumpy.encoded_array import EncodedArray
from bionumpy.encodings import BaseEncoding
import numpy as np


def chunk_from_text(text):
    return EncodedArray(np.frombuffer(bytes(text, encoding="utf8"), dtype=np.uint8), BaseEncoding)


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
"""    , "fasta": """\
>header
CTTGTTGA
>header2
CGG
"""
    , "multiline_fasta":"""\
>header
CTTGCC
GCCTCC
>header2
CCCCCC
GGGCCC
TTT
"""
    , "bed": """\
chr1\t1\t3\t.\t0\t-
chr1\t40\t60\t.\t1\t+
chr20\t400\t600\t.\t2\t+
"""
    , "vcf": """\
chr1	88362	rs4970378	A	G	.	.	.
chr1	887560	rs3748595	A	C	.	.	.
chr2	8878	rs3828047	A	G	.	.	.
"""
    , "vcf2": """\
chr1	88362	rs4970378	A	G	.	.	.
chr1	887560	rs3748595	A	CAA	.	.	.
chr2	8878	rs3828047	AGG	C	.	.	.
"""
    , "vcf_matrix": """\
chr1	883625	rs4970378	A	G\t.\t.\t.\t.\t1|1:0,4:4:6:70,6,0	1|1:0,19:19:36:358,36,0	1|1:0,3:3:6:67,6,0	1|1:0,1:1:3:34,3,0
chr1	887560	rs3748595	A	C\t.\t.\t.\t.\t0/0:7,0:7:15:0,15,163	1/1:0,30:30:81:888,81,0	1/1:0,2:2:6:68,6,0	1/1:0,1:1:3:36,3,0
chr1	887801	rs3828047	A	G\t.\t.\t.\t.\t./.	1/1:0,17:17:39:398,39,0	1/1:0,3:3:9:102,9,0	1/1:0,1:1:3:34,3,0
"""
    , "gfa_sequence": """\
S\tid1\tAACCTTGG
S\tid4\tACTG
""", "gff_entry": """\
CHROMOSOME_I	Allele	substitution	10017380	10017380	.	+.	aachange=A to T;consequence=Missense;interpolated_map_position=4.49151;public_name=q504;substitution=G/A;variation=WBVar00241143
CHROMOSOME_I	Allele	substitution	10573196	10573196	.	+.	aachange=A to T;consequence=Missense;interpolated_map_position=5.05579;public_name=vc56;substitution=G/A;variation=WBVar00275020
CHROMOSOME_I	Allele	substitution	13447169	13447169	.	+.	aachange=A to T;consequence=Missense;interpolated_map_position=20.7564;public_name=ar197;strain=GS3412;substitution=G/A;variation=WBVar00000227
CHROMOSOME_I	Allele	substitution	13721389	13721389	.	+.	aachange=A to V;consequence=Missense;interpolated_map_position=21.6433;public_name=mg334;substitution=G/A;variation=WBVar00088937
""",
    "sam": """\
@HD	VN:1.0	SO:unsorted
@SQ	SN:test_ref	LN:17637
SRR1524970.144283	16	test_ref	1706	255	25M	*	0	0	TGCTGATGAAGCAGAACAACTTTAA	]YG[^baaaa^W`ab]]````aaba	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:25	YT:Z:UU
SRR1524970.316478	16	test_ref	1706	255	24M	*	0	0	TGCTGATGAAGCAGAACAACTTTA	`\X_`aaaaaY]``b_aa_aaaaa	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:24	YT:Z:UU
""",
    "bed12": """\
chr21 10079666  10120808   uc002yiv.1  0  -  10081686  10120608  0     4   528,91,101,215, 0,1930,39750,40927,
chr21 10080031  10081687   uc002yiw.1  0  -  10080031  100800310\t0     2   200,91,    0,1565,
"""
}

buffers = {key: chunk_from_text(val) for key, val in buffer_texts.items()}

data = {
    "bed": [
        Bed6.single_entry("chr1", 1, 3, ".", "0", "-"),
        Bed6.single_entry("chr1", 40, 60, ".", "1", "+"),
        Bed6.single_entry("chr20",  400, 600, ".", "2", "+")],
    "vcf2": [
        VCFEntry.single_entry("chr1",	88361, "rs4970378",	"A",	"G", ".", ".", "."),
        VCFEntry.single_entry("chr1",	887559, "rs3748595",	"A",	"CAA", ".", ".", "."),
        VCFEntry.single_entry("chr2",	8877, "rs3828047",	"AGG",	"C", ".", ".", ".")],
    "vcf": [
        VCFEntry.single_entry("chr1",	88361, "rs4970378",	"A",	"G", ".", ".", "."),
        VCFEntry.single_entry("chr1",	887559, "rs3748595",	"A",	"C", ".", ".", "."),
        VCFEntry.single_entry("chr2",	8877, "rs3828047",	"A",	"G", ".", ".", ".")],
    "fastq": [
        SequenceEntryWithQuality.single_entry("headerishere", "CTTGTTGA", "".join("!" for _ in "CTTGTTGA")),
        SequenceEntryWithQuality.single_entry("anotherheader", "CGG", "".join("~" for _ in "CGG"))],
    "fasta": [
        SequenceEntry.single_entry("header", "CTTGTTGA"),
        SequenceEntry.single_entry("header2", "CGG")],
    "multiline_fasta": [
        SequenceEntry.single_entry("header", "CTTGCCGCCTCC"),
        SequenceEntry.single_entry("header2", "CCCCCCGGGCCCTTT")],
    "gfa_sequence": [
        SequenceEntry.single_entry("id1", "AACCTTGG"),
        SequenceEntry.single_entry("id4", "ACTG")],
    "sam": [
        SAMEntry.single_entry("SRR1524970.144283", 16, "test_ref", 1705, 255, "25M",	"*", 0, 0, "TGCTGATGAAGCAGAACAACTTTAA", "]YG[^baaaa^W`ab]]````aaba"),
        SAMEntry.single_entry("SRR1524970.316478", 16, "test_ref", 1705, 255, "24M", "*", 0, 0, "TGCTGATGAAGCAGAACAACTTTA", 	"`\X_`aaaaaY]``b_aa_aaaaa")],
    "bed12": [
        Bed12.single_entry("chr21", 10079666,  10120808,   "uc002yiv.1", "0", "-", 10081686, 10120608,  "0",     4,   [[528,91,101,215]], [[0,1930,39750,40927]]),
        Bed12.single_entry("chr21", 10080031,  10081687,   "uc002yiw.1",  "0",  "-",  10080031,  10080031,  "0",     2,   [[200,91]],    [[0,1565]])]
}


buffer_type = {"bed": Bed6Buffer,
               "vcf2": VCFBuffer,
               "vcf": VCFBuffer,
               "fastq": FastQBuffer,
               "fasta": TwoLineFastaBuffer,
               "gfa_sequence": GfaSequenceBuffer,
               "multiline_fasta": MultiLineFastaBuffer,
               "bed12": Bed12Buffer}


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
chr1	88362	rs4970378	A	G	.	.	.
chr1	887560	rs3748595	A	C	.	.	.
chr2	8878	rs3828047	A	G	.	.	.
"""
    return chunk_from_text(t)

@pytest.fixture
def vcf_buffer2():
    t = """\
chr1	88362	rs4970378	A	G	.	.	.
chr1	887560	rs3748595	A	CAA	.	.	.
chr2	8878	rs3828047	AGG	C	.	.	.
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
S\tid1\tAACCTTGG
S\tid4\tACTG
"""
    return chunk_from_text(t)

big_fastq_text = """\
@2fa9ee19-5c51-4281-abdd-eac8663f9b49 runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=109831 ch=33 start_time=2019-09-12T12:05:03Z
CGGTAGCCAGCTGCGTTCAGTATGGAAGATTTGATTTGTTTAGCGATCGCCATACTACCGTGACAAGAAAGTTGTCAGTCTTTGTGACTTGCCTGTCGCTCTATCTTCCAGACTCCTTGGTCCGTGTTCAATCCCGGTAGTAGCGACGGGCGGTGTATGTATTATCAGCGCAACAGAAACAAAGACACC
+
+&&-&%$%%$$$#)33&0$&%$''*''%$#%$%#+-5/---*&&%$%&())(&$#&,'))5769*+..*&(+28./#&1228956:7674';:;80.8>;91;>?B=%.**==?(/'($$$$*'&'**%&/));807;3A=;88>=?9498++0%"%%%%'#&5/($0.$2%&0'))*'%**&)(.%&&
@1f9ca490-2f25-484a-8972-d60922b41f6f runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=106343 ch=28 start_time=2019-09-12T12:05:07Z
GATGCATACTTCGTTCGATTTCGTTTCAACTGGACAACCTACCGTGACAAAGAAAGTTGTCGATGCTTTGTGACTTGCTGTCCTCTATCTTCAGACTCCTTGGTCCATTTCAAGACCAAACAATCAGTAGTAGCGACGGGCGGTGTGGCAATATCGCTTTCAACGAAACACAAAGAAT
+
&%&%''&'+,005<./'%*-)%(#$'$#$%&&'('$$..74483=.0412$*/,)9/194/+('%%(+1+'')+,-&,>;%%.*@@D>>?)3%%296070717%%16;<'<236800%(,734(0$7769879@;?8)09:+/4'1+**7<<4.4,%%(.)##%&'(&&%*++'&#%$
@06936a64-6c08-40e9-8a10-0fbc74812c89 runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=83531 ch=23 start_time=2019-09-12T12:03:50Z
GTTTTGTCGCTGCGTTCAGTTTATGGGTGCGGGTGTTATGATGCTTCGCTTTACGTGACAAGAAAGTTAGTAGATTGTCTTTATGTTTCTGTGGTGCTGATATTGCCACACCGCCCGATAGCTCTACCGATTGAAACACGGACCAAGGAATCGGAAATGTAGGCGAGCAGGCCGTCCTGAACACCCATTAACTTTCTTGTC
+
$&'((&%$$$.$2/=-*#'.2'&&##$$#$#&&(&+-%'(%&#"###""$$%#)%,+)+&'(,&%*((%%&%$%'+),,+,,&%$')1+*$.&+6*+(*%(&'*(''&%*+,*)('%#$$$%,$&&'&)))12)*&/*,364$%$%))$'')#%%&%$#$%$$#('$(%$%$%%$$*$&$%)''%%$$&'&$)+2++,)&%
@d6a555a1-d8dd-4e55-936f-ade7c78d9d38 runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=112978 ch=97 start_time=2019-09-12T12:03:49Z
CGTATGCTTTGAGATTCATTCAGGAGGCGGGTATTTGCTCGATCATACCATACGTGGCAAGAAAGTTGTCAGTGTCTTTGTGTTTCTCTGTGGTGCGCGATATTGCCACGCCCGTCGCTACACCGATTGAAACACGGACCGAAGTCTGAAGATAGAGCGACGAGCGAAGTCACAAAGGAACTAGAGCAACTTTTTATC
+
#$%%%%''(($$%$*-&%$%)%*'%(+($)(%$,.)##$&$$#$$&('(%&%%%%#$$%(&*('('+18/(6?65510+))'--*&&$$$,*+;/+%%&&''13&%&%(133<;9=/.2*$*657,0*&(237'85;A1/$$%'7:;;:<2:..%$)0,*.)(1)1&&1+-$$,-&(-&&####%%98:AHFEB4(%,
"""
