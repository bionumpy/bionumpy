import pytest
import bionumpy as bnp
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

@pytest.fixture()
def sam_text():
    return text.replace(' ', '\t')

@pytest.fixture
def tmp_path():
    from pathlib import Path
    path = Path('tmp_folder')
    path.mkdir(exist_ok=True)
    return path

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
    print(d)
    print(d.flag.dtype)
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

