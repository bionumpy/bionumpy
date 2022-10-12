import numpy as np
import bionumpy as bnp
from bionumpy.mutation_signature import count_mutation_types
from bionumpy.util import get_snps
from bionumpy.groupby import groupby
import logging
logging.basicConfig(level="INFO")


def main(vcf_filename, fasta_filename, flank=1):
    variants = bnp.open(vcf_filename).read_chunks()
    snps = get_snps(variants)
    snps = groupby(snps, "chromosome")
    reference = bnp.open(fasta_filename)
    counts = count_mutation_types(snps, reference, flank)
    print(counts)


def test():
    main("example_data/few_variants.vcf", "example_data/small_genome.fa.fai")


if __name__ == "__main__":
    test()
