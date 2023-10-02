"""
Example with jaccard all vs all on bedfiles.
Reads all bedfiles into memory.
Files should already be sorted.
"""

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
            print(file1, file2)
            if file1 == file2:
                continue

            run_id = frozenset((file1, file2))
            if run_id in results:
                continue

            results[run_id] = jaccard(global_size, data1, data2)

    # todo: Write results to a file
    return results


if __name__ == "__main__":
    chrom_sizes = sys.argv[1]
    files = sys.argv[2:]
    # chrom_sizes = "example_data/hg38_unix_sorted.chrom.sizes"
    # files = "benchmarks/results/bed_files/ENCFF120PGJ.sorted.bed  benchmarks/results/bed_files/ENCFF405ZPR.sorted.bed  benchmarks/results/bed_files/ENCFF539MIO.sorted.bed  benchmarks/results/bed_files/ENCFF965PER.sorted.bed benchmarks/results/bed_files/ENCFF193LLN.sorted.bed  benchmarks/results/bed_files/ENCFF410SWS.sorted.bed".split()
    result = jaccard_all_vs_all(chrom_sizes, files)
    for r in result:
        print(r)

# chrom_sizes = "example_data/hg38_unix_sorted.chrom.sizes"
# #files = "benchmarks/results/bed_files/ENCFF120PGJ.sorted.bed  benchmarks/results/bed_files/ENCFF405ZPR.sorted.bed  benchmarks/results/bed_files/ENCFF539MIO.sorted.bed  benchmarks/results/bed_files/ENCFF965PER.sorted.bed benchmarks/results/bed_files/ENCFF193LLN.sorted.bed  benchmarks/results/bed_files/ENCFF410SWS.sorted.bed".split()
# 
# files = "benchmarks/results/bed_files/ENCFF120PGJ.sorted.bed  benchmarks/results/bed_files/ENCFF405ZPR.sorted.bed  benchmarks/results/bed_files/ENCFF539MIO.sorted.bed  benchmarks/results/bed_files/ENCFF965PER.sorted.bed benchmarks/results/bed_files/ENCFF193LLN.sorted.bed  benchmarks/results/bed_files/ENCFF410SWS.sorted.bed  benchmarks/results/bed_files/ENCFF625QHR.sorted.bed benchmarks/results/bed_files/ENCFF328GQG.sorted.bed  benchmarks/results/bed_files/ENCFF522HZT.sorted.bed  benchmarks/results/bed_files/ENCFF816AEF.sorted.bed".split()
# 
# 
# jaccard_all_vs_all(chrom_sizes, files)
