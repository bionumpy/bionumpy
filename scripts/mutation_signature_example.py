import bionumpy as bnp
import typer
from bionumpy.streams.multistream import MultiStream
from bionumpy.util.cli import run_as_commandline
from bionumpy.variants.mutation_signature import count_mutation_types_genomic, count_mutation_types
from bionumpy.io.delimited_buffers import PhasedVCFMatrixBuffer
from bionumpy.io.matrix_dump import matrix_to_csv
from bionumpy import Genome
import logging
logging.basicConfig(level="INFO")


def main(vcf_filename: str, fasta_filename: str, out_filename: str = None, flank: int = 1):
    genome = Genome.from_file(fasta_filename)
    variants = genome.read_locations(vcf_filename)
    counts = count_mutation_types_genomic(variants, genome.read_sequence())
    output = matrix_to_csv(counts.counts, header=counts.alphabet)
    if out_filename is not None:
        open(out_filename, "wb").write(bytes(output.raw()))
    else:
        print(output.to_string())


def main_old(vcf_filename: str, fasta_filename: str, out_filename: str = None, flank: int = 1, genotyped: bool = False):
    if genotyped:
        variants = bnp.open(vcf_filename, buffer_type=PhasedVCFMatrixBuffer).read_chunks()
    else:
        variants = bnp.open(vcf_filename).read_chunks()

    reference = bnp.open_indexed(fasta_filename)
    sequence_lengths = reference.get_contig_lengths()
    multistream = MultiStream(sequence_lengths,
                              variants=variants,
                              reference=reference)
    counts = count_mutation_types(multistream.variants, multistream.reference, flank)
    output = matrix_to_csv(counts.counts, header=counts.alphabet)
    if out_filename is not None:
        open(out_filename, "wb").write(bytes(output.raw()))
    else:
        print(output.to_string())


def test():
    main("example_data/few_variants.vcf", "example_data/small_genome.fa")


if __name__ == "__main__":
    typer.run(main)
    # run_as_commandline(main)
