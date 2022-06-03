import numpy as np
import bionumpy as bnp
from bionumpy.delimited_buffers import VCFMatrixBuffer
from bionumpy.mutation_signature import (
    get_kmers, MutationSignatureEncoding)
from bionumpy.util import get_snps, filter_on_intervals
from bionumpy.intervals import sort_intervals, merge_intervals

import logging
logging.basicConfig(level="INFO")


def print_file(counts, flank):
    print(",".join(MutationSignatureEncoding(2*flank+1).to_string(i)
                   for i in range(counts.shape[-1])))
    for row in np.atleast_2d(counts):
        print(",".join(str(c) for c in row))


def simple_main(vcf_filename, fasta_filename, flank, do_matrix=False, bed_filename=None):
    if do_matrix:
        variants = bnp.open(vcf_filename, buffer_type=VCFMatrixBuffer)
    else:
        variants = bnp.open(vcf_filename)
    snps = get_snps(variants)
    if bed_filename:
        intervals = bnp.open(bed_filename, remove_chr=True)
        intervals = merge_intervals(sort_intervals(intervals))
        snps = filter_on_intervals(snps, intervals)
    reference = bnp.open(fasta_filename, remove_chr=True)
    counts = get_kmers(snps, reference, flank)
    print_file(counts, flank)


if __name__ == "__main__":
    import sys
    vcf_filename, fasta_filename, flank = sys.argv[1:4]
    do_matrix = len(sys.argv) > 4 and int(sys.argv[4]) == 1
    bed_filename = None if len(sys.argv) <= 5 else sys.argv[5]
    simple_main(vcf_filename, fasta_filename, int(flank),
                do_matrix, bed_filename)
