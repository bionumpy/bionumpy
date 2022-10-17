import numpy as np
import bionumpy as bnp
from bionumpy.cli import run_as_commandline
from bionumpy.mutation_signature import count_mutation_types
from bionumpy.groupby import groupby
from bionumpy.delimited_buffers import VCFMatrixBuffer
from bionumpy.io.matrix_dump import matrix_to_csv
import logging
logging.basicConfig(level="INFO")


def main(vcf_filename: str, fasta_filename: str, out_filename: str = None, flank: int = 1, genotyped: bool = True):
    if genotyped:
        variants = bnp.open(vcf_filename, buffer_type=VCFMatrixBuffer).read_chunks()
    else:
        variants = bnp.open(vcf_filename).read_chunks()
    variants = groupby(variants, "chromosome")
    reference = bnp.open(fasta_filename)
    counts = count_mutation_types(variants, reference, flank)
    output = matrix_to_csv(counts.counts, header=counts.alphabet)# , row_names=counts.row_names)
    if out_filename is not None:
        open(out_filename, "wb").write(bytes(output.raw()))
    else:
        print(output.to_string())


def test():
    main("example_data/genotype_variants.vcf", "example_data/small_genome.fa.fai")


if __name__ == "__main__":
    run_as_commandline(main)
