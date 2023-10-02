import sys
import numpy as np
import bionumpy as bnp

# Filter a vcf file by min allele count


def filter_on_allele_count(chunk, min_ac=10):
    # we work direclty with the raw encoding
    # we know that the phased encoding is 0 for 0|0, 1 for 0|1, 2 for 1|0 and 3 for 1|1
    genotypes = chunk.genotypes.raw()
    allele_counts = np.sum((genotypes == 1) + (genotypes == 2) + 2 * (genotypes == 3), axis=-1)
    return chunk[allele_counts >= min_ac]


def filter_file_on_allele_count(input_file, output_file, min_ac=10):
    f = bnp.open(input_file, buffer_type=bnp.io.PhasedVCFMatrixBuffer)
    output_file = bnp.open(output_file, "w", buffer_type=bnp.PhasedVCFMatrixBuffer)
    chunks = f.read_chunks(min_chunk_size=200000000)  # [f.read()]
    for chunk in chunks:
        print("chunk")
        filtered = filter_on_allele_count(chunk, min_ac=min_ac)
        output_file.write(filtered)


def test():
    filter_file_on_allele_count("example_data/variants_phased.vcf", "test.vcf.tmp", min_ac=1)


if __name__ == "__main__":
    filter_file_on_allele_count(sys.argv[1], sys.argv[2], int(sys.argv[3]))
