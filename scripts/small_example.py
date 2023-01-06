"""
This file contains small examples for showcasing BioNumPy
"""
import bionumpy as bnp
import numpy as np
import plotly.express as px

from bionumpy.sequence import PWM, count_kmers
from bionumpy.sequence.position_weight_matrix import PositionWeightMatrix, _pwm_from_counts



def sequence_matching():
    reads = bnp.open("example_data/reads.fq.gz").read()
    matches = bnp.sequence.match_string(reads.sequence, "AC")
    # matches is a boolean RaggedArray, so we can sum over the rows to get
    # number of matches for each read, e.g:
    # matches_per_read = np.sum(matches, axis=1)

    # .. or take a mean over the columns to see where the matches typically are
    matches_per_base = np.mean(matches, axis=0)
    fig = px.line(matches_per_base[0:150])  # plot first 150 bases
    # fig.show()


def find_kmers_in_reads():
    f = bnp.open("example_data/reads.fq.gz")
    kmer_counts = sum(
        count_kmers(chunk.sequence, k=5)
        for chunk in f.read_chunks())
    most_common = kmer_counts.most_common(10)
    return px.bar(x=most_common.labels,
                  y=most_common.counts)


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
    from bionumpy.io.motifs import read_motif
    # Read the alphabet and counts from jaspar file
    pwm = read_motif("example_data/MA0080.1.jaspar")
    reads = bnp.open("example_data/reads.fq.gz").read()
    scores = bnp.get_motif_scores(reads.sequence, pwm)
    max_scores = scores.max(axis=-1)
    return px.histogram(max_scores, nbins=15)


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
    motif_matching()
    forbes()
