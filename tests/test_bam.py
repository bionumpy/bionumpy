import numpy as np
from numpy.testing import assert_array_equal

import bionumpy as bnp
import pytest
from bionumpy.io.bam import BamIntervalBuffer
from bionumpy.util.testing import assert_encoded_raggedarray_equal
from tests.util import get_file_name



def test_read_acceptance(data_path):
    filename = data_path / "test.bam"
    f = bnp.open(filename)
    d = f.read()
    print(d)
    print(d.flag.dtype)


def test_read_intervals_acceptance(data_path):
    filename = data_path / "test.bam"
    f = bnp.open(filename, buffer_type=BamIntervalBuffer)
    d = f.read()
    print(d.start)
    assert d.start[0] == 7512371
    print(d)

@pytest.fixture()
def bam_entries(data_path):
    filename = data_path / 'small_alignments.bam'
    entries = bnp.open(filename).read()
    return entries

#@pytest.mark.xfail
def test_read_bam(bam_entries):
    start = bam_entries[:4]
    assert_array_equal(start.position, [523205, 3837782, 907877, 260353])
    assert_encoded_raggedarray_equal(start.chromosome, ['contig28', 'contig14', 'contig23', 'contig11'])
    assert_encoded_raggedarray_equal(start.sequence[:2], ['CNATCTCTTTCGTACGAGTATTTCGCGTTCTTGAGGTGAGCCTGTTAAGATCCAAATCGTTAAATAGCCGATTTCGGCTCTCGCAGTAAATTTTATAGCCATCACCTTTTCATCAATCAGCTCGCACGGCTCTACGAACCTTCGAGTTCAC',
                                                          'TTTGGCGTTAGCCACGTTTCTGACGTATAAAATGAAGCCGAGAAATCGAATCGCTGATTGCTTCATGCATCTATCATATGCCGCTGAAGAACGAGGGATCGTATGCAGCTTTTACTTTCTCAAGAACGAACGTCGGCTATTGGCTGTTTTA'])
    assert_encoded_raggedarray_equal(start.name[:2], ['ERR6054981.1', 'ERR6054981.1'])


def test_index_bam(bam_entries):
    mask = bam_entries.mapq == 60
    filtered = bam_entries[mask]
    assert_array_equal(filtered.position[:4], [523205, 3837782, 907877, 406696])


def test_write_bam(bam_entries, tmp_path):
    subset = bam_entries[bam_entries.mapq == 60]

    output_file = tmp_path / 'tmp.bam'
    with bnp.open(output_file, mode='w') as f:
        f.write(subset)
    assert open(output_file, 'rb').read()[-28:] == b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'
    new_entries = bnp.open(output_file).read()
    assert_array_equal(new_entries.position, subset.position)


def test_write_bam_after_change(bam_entries):
    bam_entries = bnp.replace(bam_entries, position=np.zeros_like(bam_entries.position))
    print(bam_entries)
    with bnp.open('tmp.bam', mode='w') as f:
        with pytest.raises(ValueError):
            f.write(bam_entries)

