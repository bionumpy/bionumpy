"""
Example with jaccard all vs all on bedfiles.
Reads all bedfiles into memory.
Files should already be sorted.
"""

import bionumpy as bnp
from bionumpy.arithmetics.similarity_measures import jaccard
from bionumpy.arithmetics import intersect
import numpy as np
import sys


def jaccard_all_vs_all(chrom_sizes, bed_file_names):
    chrom_sizes = bnp.open(chrom_sizes).read()
    all_data = {file_name: bnp.open(file_name).read() for file_name in bed_file_names}
    print(len(all_data), " files")
    results = {}

    for file1, data1 in all_data.items():
        for file2, data2 in all_data.items():
            print(file1, file2)
            if file1 == file2:
                continue

            run_id = frozenset((file1, file2))
            if run_id in results:
                continue

            #results[run_id] = intersect(data1, data2)
            results[run_id] = jaccard(chrom_sizes, data1, data2)

    print(results)
    print("Found %d pairs" % len(results))
    # todo: Write results to a file
    return results


if __name__ == "__main__":
    chrom_sizes = sys.argv[1]
    files = sys.argv[2:]
    jaccard_all_vs_all(chrom_sizes, files)


