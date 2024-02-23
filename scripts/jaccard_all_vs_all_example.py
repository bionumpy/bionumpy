"""
Calculates the jaccard similarity between all pairs of bed files.
Writes the results to stdout.
"""
import itertools
from typing import List, Tuple, Dict

import bionumpy as bnp
from bionumpy.arithmetics.similarity_measures import jaccard
from bionumpy.arithmetics import intersect
from bionumpy.genomic_data.global_offset import GlobalOffset
import numpy as np
import sys


def jaccard_func(mask_a, mask_b):
    '''Jaccard = intersection / union'''
    return (mask_a & mask_b).sum() / (mask_a | mask_b).sum()


def jaccard(chrom_sizes_file: str, bed_files: List[str]) -> Dict[Tuple[str, str], float]:
    genome = bnp.Genome.from_file(chrom_sizes_file)
    masks = {filename: genome.read_intervals(filename).get_mask() for filename in bed_files}
    results = {(file_a, file_b): jaccard_func(masks[file_a], masks[file_b])
               for file_a, file_b in itertools.combinations(masks.keys(), r=2)}
    return results


if __name__ == "__main__":
    chrom_sizes = sys.argv[1]
    files = sys.argv[2:]
    result = jaccard(chrom_sizes, files)
    for r in result.items():
        print(r)
