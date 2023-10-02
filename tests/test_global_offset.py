import pytest
from bionumpy.datatypes import Interval, ChromosomeSize
from bionumpy.genomic_data.global_offset import GlobalOffset, global_encoding
from bionumpy.util.testing import assert_bnpdataclass_equal


@pytest.fixture
def chrom_sizes():
    return ChromosomeSize.from_entry_tuples([
        ("chr1", 10),
        ("chr2", 100),
        ("chr3", 1000)])


@pytest.fixture
def local_interval():
    return Interval.from_entry_tuples([
        ("chr1", 2, 5),
        ("chr3", 15, 20),
        ("chr2", 3, 4)])


@pytest.fixture
def global_interval():
    i = Interval.from_entry_tuples([
        ("global", 2, 5),
        ("global", 125, 130),
        ("global", 13, 14)])
    i.chromosome = global_encoding.encode(i.chromosome)
    return i


def test_global_offset(local_interval, global_interval, chrom_sizes):
    global_offset = GlobalOffset(chrom_sizes)
    new_interval = global_offset.from_local_interval(local_interval)
    print(new_interval)
    assert_bnpdataclass_equal(new_interval, global_interval)
