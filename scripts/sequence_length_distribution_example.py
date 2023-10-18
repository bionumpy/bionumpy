import gzip

from numpy.testing import assert_array_equal

import bionumpy as bnp
import numpy as np
import sys

import bionumpy.io.fastq_buffer
from bionumpy.io.parser import NumpyFileReader


def get_sorted_sequence_length_distribution(file_name):
    counts = bnp.bincount(chunk.sequence.lengths for chunk in bnp.open(file_name).read_chunks())
    sizes = np.flatnonzero(counts)
    return sizes, counts[sizes]


def test():
    get_sorted_sequence_length_distribution("example_data/big.fq.gz")

def _test_2():
    sizes, counts = get_sorted_sequence_length_distribution("../benchmarks/results/dna_sequences/ENCFF689IPX.fq.gz")
    assert_array_equal(sizes, [32])
    assert_array_equal(counts, [8886714])

if __name__ == "__main__":
    file_name = sys.argv[1]
    out_file = sys.argv[2]

    sizes, counts = get_sorted_sequence_length_distribution(file_name)
    with open(out_file, "w") as f:
        f.write("\n".join(["%d %d" % (size, count) for size, count in zip(sizes, counts)]) + "\n")
