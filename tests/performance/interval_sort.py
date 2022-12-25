import pytest
import dataclasses
import bionumpy as bnp
import numpy as np
from bionumpy.util.testing import assert_bnpdataclass_equal
from bionumpy.encodings.string_encodings import StringEncoding

@pytest.fixture
def intervals():
    return bnp.open("example_data/ctcf.bed.gz").read()

@pytest.fixture
def chrom_sizes():
    return bnp.open("example_data/hg38.chrom.sizes").read()

def sort_func(intervals):
    args = np.lexsort((intervals.start, intervals.chromosome))
    return intervals[args]


def _test_interval_sort(intervals, chrom_sizes):
    truth = bnp.arithmetics.sort_intervals(intervals, sort_order=chrom_sizes.name.tolist())
    encoding = StringEncoding(chrom_sizes.name)
    new_intervals = dataclasses.replace(intervals,
                                        chromosome=encoding.encode(intervals.chromosome))
    result = sort_func(new_intervals)
    np.testing.assert_array_equal(truth.start, result.start)


def test_interval_sort_benchmark(benchmark, intervals, chrom_sizes):
    truth = bnp.arithmetics.sort_intervals(intervals, sort_order=chrom_sizes.name.tolist())
    encoding = StringEncoding(chrom_sizes.name)
    new_intervals = dataclasses.replace(intervals,
                                        chromosome=encoding.encode(intervals.chromosome))
    result = benchmark(sort_func, new_intervals)
    np.testing.assert_array_equal(truth.start, result.start)
    
    
    
