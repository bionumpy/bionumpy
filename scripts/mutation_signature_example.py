import numpy as np
import bionumpy as bnp
from bionumpy.streams.multistream import MultiStream
from bionumpy.cli import run_as_commandline
from bionumpy.mutation_signature import count_mutation_types
from bionumpy.io.delimited_buffers import PhasedVCFMatrixBuffer
from bionumpy.io.matrix_dump import matrix_to_csv
import logging
logging.basicConfig(level="INFO")


def main(vcf_filename: str, fasta_filename: str, out_filename: str = None, flank: int = 1, genotyped: bool = False):
    if genotyped:
        variants = bnp.open(vcf_filename, buffer_type=PhasedVCFMatrixBuffer).read_chunks()
    else:
        variants = bnp.open(vcf_filename).read_chunks()

    reference = bnp.open_indexed(fasta_filename)
    sequence_lengths = reference.get_contig_lengths()
    #sequence_lengths = {name: sequence_lengths[name] for name in 
    #                     sorted(sequence_lengths.keys(), key=alpha_numeric_key_func)}
    multistream = MultiStream(sequence_lengths,
                              variants=variants,
                              reference=reference)
    # multistream.set_key_functions(variants=lambda x: x if x.startswith("chr") else "chr"+x)
    counts = count_mutation_types(multistream.variants, multistream.reference, flank)
    output = matrix_to_csv(counts.counts, header=counts.alphabet)# , row_names=counts.row_names)
    if out_filename is not None:
        open(out_filename, "wb").write(bytes(output.raw()))
    else:
        print(output.to_string())


def test():
    main("example_data/few_variants.vcf", "example_data/small_genome.fa")


if __name__ == "__main__":
    run_as_commandline(main)
