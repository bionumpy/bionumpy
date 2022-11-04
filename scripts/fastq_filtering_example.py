import numpy as np
import bionumpy as bnp
from bionumpy import streamable


@streamable()
def filter_reads_on_mean_base_quality(reads, minimum_base_quality=20):
    mask = np.mean(reads.quality, axis=-1) > minimum_base_quality
    return reads[mask]


@streamable()
def filter_reads_on_minimum_base_quality(reads, min_base_quality=5):
    mask = np.min(reads.quality, axis=-1) > min_base_quality
    return reads[mask]


def test(file="example_data/big.fq.gz"):
    reads = bnp.open(file).read_chunks()
    reads = filter_reads_on_mean_base_quality(reads, 10)
    reads = filter_reads_on_minimum_base_quality(reads, 1)
    print("Number of reads after filtering: ", streamable(sum)(len)(reads))


if __name__ == "__main__":
    test()