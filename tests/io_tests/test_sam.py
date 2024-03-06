import pytest
import bionumpy as bnp
from bionumpy.io.headers import SamHeader
from bionumpy.util.testing import assert_encoded_array_equal, assert_bnpdataclass_equal

text = '''\
@HD VN:1.6 SO:coordinate
@SQ SN:ref LN:47
ref 516 ref 1 0 14M2D31M * 0 0 AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT *
r001 99 ref 7 30 14M1D3M = 39 41 TTAGATAAAGGATACTG *
* 768 ref 8 30 1M * 0 0 * * CT:Z:.;Warning;Note=Ref wrong?
r002 0 ref 9 30 3S6M1D5M * 0 0 AAAAGATAAGGATA * PT:Z:1;4;+;homopolymer
r003 0 ref 9 30 5H6M * 0 0 AGCTAA * NM:i:1
r004 0 ref 18 30 6M14N5M * 0 0 ATAGCTTCAGC *
r003 2064 ref 31 30 6H5M * 0 0 TAGGC * NM:i:0
r001 147 ref 39 30 9M = 7 -41 CAGCGGCAT * NM:i:1
'''

header_text = '''\
@HD	VN:1.6	SO:queryname
@SQ	SN:contig0	LN:4078262
@SQ	SN:contig1	LN:4078262
@SQ	SN:contig2	LN:4078262
@SQ	SN:contig3	LN:4078262
@SQ	SN:contig4	LN:4077762
@SQ	SN:contig5	LN:4078262
@SQ	SN:contig6	LN:4078262
@SQ	SN:contig7	LN:4078263
@SQ	SN:contig8	LN:3914785
@SQ	SN:contig9	LN:3914785
@SQ	SN:contig10	LN:3914785
@SQ	SN:contig11	LN:3914785
@SQ	SN:contig12	LN:3913785
@SQ	SN:contig13	LN:3914285
@SQ	SN:contig14	LN:3913785
@SQ	SN:contig15	LN:3914791
@SQ	SN:contig16	LN:2552982
@SQ	SN:contig17	LN:2552482
@SQ	SN:contig18	LN:2552982
@SQ	SN:contig19	LN:2552982
@SQ	SN:contig20	LN:2552782
@SQ	SN:contig21	LN:2551326
@SQ	SN:contig22	LN:2552982
@SQ	SN:contig23	LN:2552982
@SQ	SN:contig24	LN:2552982
@SQ	SN:contig25	LN:2552989
@SQ	SN:contig26	LN:2818356
@SQ	SN:contig27	LN:2818356
@SQ	SN:contig28	LN:2818356
@SQ	SN:contig29	LN:2818156
@SQ	SN:contig30	LN:2816556
@SQ	SN:contig31	LN:2817856
@SQ	SN:contig32	LN:2818361
@SQ	SN:contig33	LN:9242430
@SQ	SN:contig34	LN:9239830
@SQ	SN:contig35	LN:2225803
@SQ	SN:contig36	LN:2225303
@SQ	SN:contig37	LN:2225603
@SQ	SN:contig38	LN:2225203
@SQ	SN:contig39	LN:2225303
@SQ	SN:contig40	LN:2225803
@SQ	SN:contig41	LN:2225803
@SQ	SN:contig42	LN:2225806
@SQ	SN:contig43	LN:8294940
@SQ	SN:contig44	LN:8294941
@SQ	SN:contig45	LN:1768679
@SQ	SN:contig46	LN:1768679
@SQ	SN:contig47	LN:1768679
@SQ	SN:contig48	LN:1768479
@SQ	SN:contig49	LN:1768482
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -t 8 -5SPM data/athalia_rosea/real/big/10/10000000/1/not_assembled/50/0/0.0/0/0.0/0.0/6000/hifiasm.hic.p_ctg.fa data/athalia_rosea/real/big/hic/10000000/1/reads1.fq.gz data/athalia_rosea/real/big/hic/10000000/1/reads2.fq.gz
@PG	ID:samtools	PN:samtools	PP:bwa	VN:1.17	CL:samtools view -buS -
@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.10	CL:samtools view -h samtools.105676.2489.tmp.0000.bam
'''

@pytest.fixture()
def sam_text():
    return text.replace(' ', '\t')


@pytest.fixture
def sam_filename(tmp_path, sam_text):
    filename = tmp_path / 'test.sam'
    filename.write_text(sam_text)
    return filename

@pytest.fixture
def sam_out_filename(tmp_path):
    filename = tmp_path / 'testout.sam'
    return filename

# @pytest.mark.xfail
def test_sam_read(sam_filename):
    f = bnp.open(sam_filename)
    d = f.read()
    assert_encoded_array_equal(d.extra[-1], 'NM:i:1')


def test_sam_write(sam_filename, sam_out_filename):
    d = bnp.open(sam_filename).read()
    with bnp.open(sam_out_filename, mode='w') as f:
        f.write(d)
    d2 = bnp.open(sam_out_filename).read()
    assert_bnpdataclass_equal(d, d2)


def test_sam_write_do(sam_filename, sam_out_filename):
    d = bnp.open(sam_filename).read().get_data_object()
    print(d)
    with bnp.open(sam_out_filename, mode='w') as f:
        f.write(d)
    d2 = bnp.open(sam_out_filename).read()
    assert_bnpdataclass_equal(d, d2)


def test_sam_header():
    header = SamHeader.from_text(header_text)
    assert header.contig_dict['contig34'] == 9239830
