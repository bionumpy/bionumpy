import numpy as np
from bionumpy import bnp_open
from bionumpy.delimited_buffers import VCFMatrixBuffer
from bionumpy.mutation_signature import get_kmers, MutationSignatureEncoding, get_snps
import logging
logging.basicConfig(level="INFO")


def print_file(counts, flank):
    print(",".join(MutationSignatureEncoding(2*flank+1).to_string(i)
                   for i in range(counts.shape[-1])))
    for row in np.atleast_2d(counts):
        print(",".join(str(c) for c in row))


def simple_main(vcf_filename, fasta_filename, flank, do_matrix=False):
    if do_matrix:
        variants = bnp_open(vcf_filename, buffer_type=VCFMatrixBuffer)
    else:
        variants = bnp_open(vcf_filename)
    snps = get_snps(variants)
    reference = bnp_open(fasta_filename, remove_chr=True)
    counts = get_kmers(snps, reference, flank)
    print_file(counts, flank)


if __name__ == "__main__":
    import sys
    vcf_filename, fasta_filename, flank = sys.argv[1:4]
    do_matrix = len(sys.argv) > 4 and int(sys.argv[4]) == 1
    simple_main(vcf_filename, fasta_filename, int(flank), do_matrix)
