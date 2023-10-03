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
# prefix = "../"
prefix = ""


def jaccard_all_vs_all(chrom_sizes, bed_file_names):
    chrom_sizes = bnp.open(prefix+chrom_sizes).read()
    global_offset = GlobalOffset(chrom_sizes)
    all_data = {file_name: global_offset.from_local_interval(bnp.open(prefix+file_name).read(), do_clip=True)
                for file_name in bed_file_names}
    global_size = {"global": sum(chrom_sizes.size)}
    results = {}

    for file1, data1 in all_data.items():
        for file2, data2 in all_data.items():
            # print(file1, file2)
            if file1 == file2:
                continue

            run_id = frozenset((file1, file2))
            if run_id in results:
                continue

            results[run_id] = jaccard(global_size, data1, data2)

    # todo: Write results to a file
    return results


def jaccard_func(mask_a, mask_b):
    return (mask_a & mask_b).sum() / (mask_a | mask_b).sum()


def jaccard(chrom_sizes_file: str, bed_files: List[str]):
    genome = bnp.Genome.from_file(chrom_sizes_file)
    masks = {filename: genome.read_intervals(filename).get_mask() for filename in bed_files}
    results = {(file_a, file_b): jaccard_func(masks[file_a], masks[file_b]) for file_a, file_b in itertools.combinations(masks.keys(), r=2)}
    return results


def _test_profiling():
    chrom, *inputs = '../example_data/hg38_unix_sorted.chrom.sizes ../benchmarks/results/intervals/ENCFF143HTO_mapped_reads_100k.bed ../benchmarks/results/intervals/ENCFF227NIG_mapped_reads_100k.bed'.split()
    jaccard(chrom, inputs)


if __name__ == "__main__":
    chrom_sizes = sys.argv[1]
    files = sys.argv[2:]
    # chrom_sizes = "example_data/hg38_unix_sorted.chrom.sizes"
    # files = "benchmarks/results/bed_files/ENCFF120PGJ.sorted.bed  benchmarks/results/bed_files/ENCFF405ZPR.sorted.bed  benchmarks/results/bed_files/ENCFF539MIO.sorted.bed  benchmarks/results/bed_files/ENCFF965PER.sorted.bed benchmarks/results/bed_files/ENCFF193LLN.sorted.bed  benchmarks/results/bed_files/ENCFF410SWS.sorted.bed".split()
    # result = jaccard_all_vs_all(chrom_sizes, files)
    result = jaccard(chrom_sizes, files)
    for r in result.items():
        print(r)

# chrom_sizes = "example_data/hg38_unix_sorted.chrom.sizes"
# #files = "benchmarks/results/bed_files/ENCFF120PGJ.sorted.bed  benchmarks/results/bed_files/ENCFF405ZPR.sorted.bed  benchmarks/results/bed_files/ENCFF539MIO.sorted.bed  benchmarks/results/bed_files/ENCFF965PER.sorted.bed benchmarks/results/bed_files/ENCFF193LLN.sorted.bed  benchmarks/results/bed_files/ENCFF410SWS.sorted.bed".split()
# 
# files = "benchmarks/results/bed_files/ENCFF120PGJ.sorted.bed  benchmarks/results/bed_files/ENCFF405ZPR.sorted.bed  benchmarks/results/bed_files/ENCFF539MIO.sorted.bed  benchmarks/results/bed_files/ENCFF965PER.sorted.bed benchmarks/results/bed_files/ENCFF193LLN.sorted.bed  benchmarks/results/bed_files/ENCFF410SWS.sorted.bed  benchmarks/results/bed_files/ENCFF625QHR.sorted.bed benchmarks/results/bed_files/ENCFF328GQG.sorted.bed  benchmarks/results/bed_files/ENCFF522HZT.sorted.bed  benchmarks/results/bed_files/ENCFF816AEF.sorted.bed".split()
# 
# 
# jaccard_all_vs_all(chrom_sizes, files)
