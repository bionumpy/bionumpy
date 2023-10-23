"""
Example with jaccard all vs all on bedfiles.
Reads all bedfiles into memory.
Files should already be sorted.
"""
import itertools
from typing import List

import bionumpy as bnp
from bionumpy.arithmetics.similarity_measures import jaccard
from bionumpy.arithmetics import intersect
from bionumpy.genomic_data.global_offset import GlobalOffset
import numpy as np
import sys


def jaccard_func(mask_a, mask_b):
    return (mask_a & mask_b).sum() / (mask_a | mask_b).sum()


def jaccard(chrom_sizes_file: str, bed_files: List[str]):
    genome = bnp.Genome.from_file(chrom_sizes_file)
    masks = {filename: genome.read_intervals(filename).get_mask() for filename in bed_files}
    results = {(file_a, file_b): jaccard_func(masks[file_a], masks[file_b])
               for file_a, file_b in itertools.combinations(masks.keys(), r=2)}
    return results


def _test_profiling():
    chrom, *inputs = '../example_data/hg38_unix_sorted.chrom.sizes ../benchmarks/results/intervals/ENCFF143HTO_mapped_reads_100k.bed ../benchmarks/results/intervals/ENCFF227NIG_mapped_reads_100k.bed'.split()
    score = list(jaccard(chrom, inputs).values())[0]
    np.testing.assert_allclose(score, 0.00884, rtol=1e-3)


if __name__ == "__main__":
    chrom_sizes = sys.argv[1]
    files = sys.argv[2:]
    result = jaccard(chrom_sizes, files)
    for r in result.items():
        print(r)
