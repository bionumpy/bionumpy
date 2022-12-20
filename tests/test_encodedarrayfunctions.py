import pytest
import dataclasses
import bionumpy as bnp
import numpy as np
from bionumpy.encodings.string_encodings import StringEncoding


@pytest.fixture
def array_1():
    return bnp.as_encoded_array("ACGTTGCA", target_encoding=bnp.DNAEncoding)


@pytest.fixture
def array_2():
    return bnp.as_encoded_array("TTTTGGGG", target_encoding=bnp.DNAEncoding)


def test_lex_sort(array_1):
    np.testing.assert_array_equal(np.lexsort((array_1,)), [0, 7, 1, 6, 2, 5, 3, 4])


def test_lex_sort_2(array_1, array_2):
    np.testing.assert_array_equal(np.lexsort((array_1, array_2)), 
                                  [7, 6, 5, 4, 0, 1, 2, 3])


@pytest.fixture
def intervals():
    return bnp.open("example_data/ctcf.bed.gz").read()

@pytest.fixture
def chrom_sizes():
    return bnp.open("example_data/hg38.chrom.sizes").read()


def test_interval_sort(intervals, chrom_sizes):
    truth = bnp.arithmetics.sort_intervals(intervals, sort_order=chrom_sizes.name.tolist())
    encoding = StringEncoding(chrom_sizes.name)
    new_intervals = dataclasses.replace(intervals,
                                        chromosome=encoding.encode(intervals.chromosome))
    result = bnp.arithmetics.sort_intervals(new_intervals)
    np.testing.assert_array_equal(truth.start, result.start)
