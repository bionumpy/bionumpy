import bionumpy as bnp
import numpy as np

import bionumpy.io.vcf_buffers


def remove_variants_with_N(chunk):
    # find and remove variants that have N in ref or alt sequence
    mask = (np.sum(chunk.ref_seq == "N", axis=1) > 0) | (np.sum(chunk.alt_seq == "N", axis=1) > 0)
    print("Found %d sequences with N in ref or var" % np.sum(mask))
    new_variants = chunk[~mask]

    # write new variants to file:
    with bnp.open("tmp_variants.vcf", "w", buffer_type=bionumpy.io.vcf_buffers.VCFMatrixBuffer) as f:
        f.write(new_variants)


def phsased_vcf_example(chunk):
    genotypes = chunk.genotypes
    print(chunk)


def test():
    analyses = [remove_variants_with_N, phsased_vcf_example]
    analyses = [remove_variants_with_N]
    for func in analyses:
        f = bnp.open("example_data/variants_with_header.vcf", buffer_type=bionumpy.io.vcf_buffers.VCFMatrixBuffer)
        data = f.read()
        func(data)


if __name__ == "__main__":
    test()
