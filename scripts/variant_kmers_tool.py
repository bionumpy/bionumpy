import bionumpy as bnp
import typer
from bionumpy.genomic_data import Genome


def main(fasta_filename: str, vcf_filename: str, k: int):
    genome = Genome.from_file(fasta_filename)
    variants = genome.read_locations(vcf_filename, has_numeric_chromosomes=True)
    windows = variants.get_windows(flank=k-1)
    sequences = genome.read_sequence()[windows]
    sequences = bnp.as_encoded_array(sequences, bnp.DNAEncoding)
    kmers = bnp.get_kmers(sequences, k)
    return kmers


if __name__ == '__main__':
    typer.run(main)
