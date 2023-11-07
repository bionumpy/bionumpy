import numpy as np
from numpy.testing import assert_array_equal

import bionumpy as bnp
import pytest
from bionumpy.io.bam import BamIntervalBuffer
from bionumpy.alignments import alignment_to_interval
from bionumpy.util.testing import assert_encoded_raggedarray_equal


#@pytest.mark.skip
def test_read_acceptance():
    filename = "example_data/test.bam"
    f = bnp.open(filename)
    d = f.read()
    print(d)
    print(d.flag.dtype)


def test_read_intervals_acceptance():
    filename = "example_data/test.bam"
    f = bnp.open(filename, buffer_type=BamIntervalBuffer)
    d = f.read()
    print(d.start)
    assert d.start[0] == 7512371
    print(d)

@pytest.fixture()
def bam_entries():
    filename = 'example_data/small_alignments.bam'
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


@pytest.mark.xfail
def test_write_bam(bam_entries):
    with bnp.open('tmp.bam', mode='w') as f:
        f.write(bam_entries[:10])
    new_entries = bnp.open('tmp.bam').read()
    assert_array_equal(new_entries.position, bam_entries[:10].position)