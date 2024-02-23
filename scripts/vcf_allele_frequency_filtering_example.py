import sys
import numpy as np
import bionumpy as bnp
from itertools import islice

# Filter a vcf file by min allele count


def filter_on_allele_count(chunk, min_ac=10):
    # we work direclty with the raw encoding
    # we know that the phased encoding is 0 for 0|0, 1 for 0|1, 2 for 1|0 and 3 for 1|1
    genotypes = chunk.genotypes.raw()
    allele_counts = np.sum((genotypes == 1) + (genotypes == 2) + 2 * (genotypes == 3), axis=-1)
    return chunk[allele_counts >= min_ac]


def filter_file_on_allele_count(input_file, output_file, min_ac=10):
    with bnp.open(output_file, "w") as output_file:
        chunks = bnp.open(input_file).read_chunks()
        for chunk in islice(chunks, 0, None):
            output_file.write(chunk[chunk.info.AC.ravel() >= min_ac])


def test():
    filter_file_on_allele_count("example_data/variants_with_header.vcf", "test.vcf", min_ac=1)


def _test_profile():
    filter_file_on_allele_count("../benchmarks/results/vcfs/big_phased.vcf.gz", "tmp.test.vcf", min_ac=10)


if __name__ == "__main__":
    filter_file_on_allele_count(sys.argv[1], sys.argv[2], int(sys.argv[3]))
