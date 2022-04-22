import numpy as np
from bionumpy import bnp_open
from bionumpy.mutation_signature import get_kmers, MutationSignatureEncoding, get_snps
import logging
logging.basicConfig(level="INFO")


def print_file(counts, flank):
    print(",".join(MutationSignatureEncoding.to_string(i, 2*flank)
                   for i in range(counts.size)))
    for row in np.atleast_2d(counts):
        print(",".join(str(c) for c in counts))


def simple_main(vcf_filename, fasta_filename, flank):
    variants = bnp_open(vcf_filename)
    reference = bnp_open(fasta_filename, remove_chr=True)
    counts = get_kmers(variants, reference, flank)
    print_file(counts, flank)
    # print('"","test"')
    # for i, count in enumerate(counts):
    #     text = MutationSignatureEncoding.to_string(i, 2*flank)
    #     print(f'"{text}",{count}')


def main(vcf_filename, bed_filename, fasta_filename):
    # bed_file = FullBedFile.from_bed_buffer_stream(
    # bedfile = BufferedNumpyParser(open(bed_filename, "rb"), BedBuffer).get_chunks()
    # intervals = SortedIntervals(next(bedfile).get_data())
    vcf = BufferedNumpyParser(open(vcf_filename, "rb"), VCFMatrixBuffer).get_chunks()
    snps, genotypes = next(vcf).get_entries()
    reference = IndexedFasta(fasta_filename, add_chr=True)["chr1"]
                                 #counts = get_kmers((snps, genotypes), None, reference, 1)
    counts = get_kmers((snps, None), None, reference, 1)
    print(",".join(MutationSignatureEncoding.to_string(c, 4) for c in range(counts.shape[-1])))
    for row in counts:
        print(",".join(str(c) for c in row))


if __name__ == "__main__":
    import sys
    # test_get_kmers()
    simple_main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
