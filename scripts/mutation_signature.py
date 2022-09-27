import numpy as np
import bionumpy as bnp
from bionumpy.mutation_signature import (
    get_kmers, MutationSignatureEncoding)
from bionumpy.util import get_snps
from bionumpy.groupby import groupby
import logging
logging.basicConfig(level="INFO")


def print_file(counts, flank):
    print(",".join(MutationSignatureEncoding(2*flank+1).to_string(i)
                   for i in range(counts.shape[-1])))
    for row in np.atleast_2d(counts):
        print(",".join(str(c) for c in row))


def main(vcf_filename, fasta_filename, flank=1):
    variants = bnp.open(vcf_filename).read_chunks()
    snps = get_snps(variants)
    snps = groupby(snps, "chromosome")
    reference = bnp.open(fasta_filename)
    counts = get_kmers(snps, reference, flank)
    print_file(counts, flank)


def test():
    main("example_data/few_variants.vcf", "example_data/small_genome.fa.fai")


if __name__ == "__main__":
    test()
