import bionumpy as bnp
from bionumpy.genomic_data import Genome


def main(fasta_filename: str, vcf_filename: str, k: int):
    # Read genome and variants
    genome = bnp.Genome.from_file(fasta_filename, filter_function=None)
    variants = genome.read_locations(vcf_filename, has_numeric_chromosomes=False)

    # get only snps from vcf
    snp_mask = (variants.get_data_field('ref_seq').shape[-1] == 1) & (variants.get_data_field('alt_seq').shape[-1] == 1)
    variants = variants[snp_mask]

    # Get windows around these snps
    windows = variants.get_windows(flank=k-1)

    # Use the windows to extra sequences (kmers) fr
    sequences = genome.read_sequence()[windows]
    sequences = bnp.as_encoded_array(sequences, bnp.DNAEncoding)
    reference_kmers = bnp.get_kmers(sequences, k)

    # Extract kmers with alternative allele on SNPs
    sequences[:, k-1] = variants.get_data_field('alt_seq').ravel()
    sequences = bnp.as_encoded_array(sequences, bnp.DNAEncoding)
    alt_kmers = bnp.get_kmers(sequences, k)

    return reference_kmers, alt_kmers


def test():
    main("example_data/sacCer3.fa", "example_data/sacCer3_sample_variants.vcf.gz", k=31)


if __name__ == '__main__':
    import typer
    typer.run(main)
