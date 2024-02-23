"""
Calculates the jaccard similarity between all pairs of bed files.
Writes the results to stdout.
"""
import itertools
from typing import List, Tuple, Dict
import bionumpy as bnp



def jaccard_func(mask_a: bnp.genomic_data.GenomicArray, mask_b: bnp.genomic_data.GenomicArray):
    '''Jaccard = intersection / union'''
    return (mask_a & mask_b).sum() / (mask_a | mask_b).sum()


def jaccard(chrom_sizes_file: str, bed_files: List[str]) -> Dict[Tuple[str, str], float]:
    genome = bnp.Genome.from_file(chrom_sizes_file)
    masks = {filename: genome.read_intervals(filename).get_mask() for filename in bed_files}
    results = {(file_a, file_b): jaccard_func(masks[file_a], masks[file_b])
               for file_a, file_b in itertools.combinations(masks.keys(), r=2)}
    return results


def test():
    chrom_sizes = 'example_data/hg38.chrom.sizes'
    files = ['example_data/ctcf.bed.gz', 'example_data/znf263.bed.gz']
    assert len(jaccard(chrom_sizes, files)) == 1


if __name__ == "__main__":
    import sys
    chrom_sizes = sys.argv[1]
    files = sys.argv[2:]
    result = jaccard(chrom_sizes, files)
    for r in result.items():
        print(r)
