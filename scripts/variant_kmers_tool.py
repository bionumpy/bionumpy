import bionumpy as bnp
import typer
from bionumpy.genomic_data import Genome


def main(fasta_filename: str, vcf_filename: str, k: int):
    genome = Genome.from_file(fasta_filename)
    variants = genome.read_locations(vcf_filename, has_numeric_chromosomes=True)
    snp_mask = (variants.get_data_field('ref_seq').shape[-1] == 1) & (variants.get_data_field('alt_seq').shape[-1] == 1)
    variants = variants[snp_mask]
    windows = variants.get_windows(flank=k-1)
    sequences = genome.read_sequence()[windows]
    sequences[:, k-1] = variants.get_data_field('alt_seq').ravel()
    sequences = bnp.as_encoded_array(sequences, bnp.DNAEncoding)
    kmers = bnp.get_kmers(sequences, k)
    return kmers


# main('/home/knut/Data/hg38.fa', '/home/knut/Data/variants.vcf.gz', 5)

if __name__ == '__main__':
    typer.run(main)
