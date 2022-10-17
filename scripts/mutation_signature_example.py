import numpy as np
import bionumpy as bnp
from bionumpy.mutation_signature import count_mutation_types
from bionumpy.groupby import groupby
from bionumpy.delimited_buffers import VCFMatrixBuffer
import logging
logging.basicConfig(level="INFO")


def main(vcf_filename, fasta_filename, flank=1, genotyped=False):
    if genotyped:
        variants = bnp.open(vcf_filename, buffer_type=VCFMatrixBuffer).read_chunks()
    else:
        variants = bnp.open(vcf_filename).read_chunks()
    variants = groupby(variants, "chromosome")
    reference = bnp.open(fasta_filename)
    counts = count_mutation_types(variants, reference, flank)
    print(counts)


def test():
    main("example_data/few_variants.vcf", "example_data/small_genome.fa.fai")


if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        main(*sys.argv[1:], genotyped=True)
    else:
        test()
    
