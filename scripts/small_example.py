"""
This file contains small examples for showcasing BioNumPy
"""
import bionumpy as bnp
import numpy as np
import plotly.express as px

from bionumpy.sequence import PWM
from bionumpy.sequence.position_weight_matrix import PositionWeightMatrix, _pwm_from_counts


def sequence_matching():
    reads = bnp.open("example_data/reads.fq.gz").read()
    matches = bnp.sequence.match_string(reads.sequence, "AC")
    # matches is a boolean RaggedArray, so we can sum over the rows to get
    # number of matches for each read
    matches_per_read = np.sum(matches, axis=1)

    # .. or take a mean over the columns to see where the matches typically are
    matches_per_base = np.mean(matches, axis=0)
    fig = px.line(matches_per_base[0:150])  # plot first 150 bases
    # fig.show()


def find_kmers_in_reads():
    f = bnp.open("example_data/reads.fq.gz")
    for chunk in f.read_chunks():
        # Change encoding for ultra-fast kmer-hashing
        sequences = bnp.change_encoding(chunk.sequence, bnp.DNAEncoding)
        kmers = bnp.get_kmers(sequences, 5)
        # kmers is now a RaggedArray that can be indexed
        first_kmer_in_each_read = kmers[:, 0]
        all_kmers_as_np_array = kmers.raw()


def filter_fastq_reads_on_base_quality():
    in_file = bnp.open("example_data/reads.fq.gz")
    out_file = bnp.open("filtered.fq.gz", "w")
    for chunk in in_file.read_chunks():
        # Create a boolean mask with reads having avg qual > 30
        keep = np.mean(chunk.quality, axis=1) > 30
        # We can index the chunk with the boolean mask
        out_file.write(chunk[keep])

        # .. or, keep reads where minimum qual > 10
        keep = np.min(chunk.quality, axis=1) > 10
        # ..


def plot_mean_base_quality_across_reads():
    reads = bnp.open("example_data/reads.fq.gz").read_chunk()
    mean_qual_per_base = bnp.mean(reads.quality, axis=0)
    fig = px.line(mean_qual_per_base[0:150])
    #fig.show()


def motif_matching():
    from bionumpy.io.jaspar import read_jaspar_matrix

    # Read the alphabet and counts from jaspar file
    alphabet, matrix = read_jaspar_matrix("example_data/MA0080.1.jaspar")

    # Convert counts to position weight matrix
    pwm = _pwm_from_counts(matrix)

    # Make an array-class for the alphabet
    arrayclass = get_alphabet_array_class(alphabet)

    # Get the motif score function
    motif_score = PositionWeightMatrix(pwm, arrayclass)

    # todo


def forbes():
    import bionumpy as bnp
    from bionumpy.arithmetics import forbes, sort_intervals

    # Read the bed files
    a = bnp.open("example_data/ctcf.bed.gz").read()
    b = bnp.open("example_data/znf263.bed.gz").read()

    # We also need chromosomes to do sorting
    chromosome_sizes = bnp.open("example_data/hg38.chrom.sizes").read()

    # sort the bed files
    a_sorted = sort_intervals(a, sort_order=chromosome_sizes.name.tolist())
    b_sorted = sort_intervals(b, sort_order=chromosome_sizes.name.tolist())

    similarity = forbes(chromosome_sizes, a_sorted, b_sorted)


def test():
    find_kmers_in_reads()
    filter_fastq_reads_on_base_quality()
    sequence_matching()
    plot_mean_base_quality_across_reads()
    #motif_matching()
    forbes()
