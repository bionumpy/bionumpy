import pytest

import bionumpy as bnp
from bionumpy.arithmetics import forbes, sort_intervals

# Read the bed files

@pytest.mark.skip
def test():
    a = bnp.open("example_data/ctcf.bed.gz").read()
    b = bnp.open("example_data/znf263.bed.gz").read()

    # We also need chromosomes to do sorting
    chromosome_sizes = bnp.open("example_data/hg38.chrom.sizes").read()
    print(chromosome_sizes)

    # sort the bed files
    a_sorted = sort_intervals(a, sort_order=chromosome_sizes.name.tolist())
    b_sorted = sort_intervals(b, sort_order=chromosome_sizes.name.tolist())

    similarity = forbes(chromosome_sizes, a_sorted, b_sorted)
    print(similarity)
