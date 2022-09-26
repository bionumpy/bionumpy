"""
This scripts performs some of the same actions that the commonly
used tools FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
"""

import numpy as np
import bionumpy as bnp
from bionumpy.npdataclassstream import streamable
import matplotlib.pyplot as plt


@streamable
def get_base_quality_histogram(reads):
    return np.bincount(reads.quality.ravel(), minlength=60)


def plot_base_qualities(reads):
    qualities = sum(get_base_quality_histogram(reads))
    plt.plot(qualities)
    plt.show()


@streamable
def get_gc_content(reads):
    sequences = reads.sequence
    mask = (sequences == ord("G")) | (sequences == ord("C"))
    return np.sum(mask, axis=-1) / sequences.shape.lengths


def plot_gc_content(reads):
    gc_content = np.concatenate(list(get_gc_content(reads)))
    plt.hist(gc_content, bins=50)
    plt.show()


@streamable
def get_quality_scores_as_matrix(reads, limit_at_n_bases=150):
    # get a padded matrix (not all reads are the same length)
    matrix = reads.quality.as_padded_matrix(side="right", fill_value=0)[:,0:limit_at_n_bases]
    return matrix


def plot_averege_quality_scores_per_base(reads):
    scores = np.mean(np.concatenate(list(get_quality_scores_as_matrix(reads)), axis=0), axis=0)
    plt.plot(scores)
    plt.show()


def test():
    examples = [plot_base_qualities, plot_gc_content, plot_averege_quality_scores_per_base]
    for example in examples:
        reads = bnp.open("example_data/big.fq.gz").read_chunks(chunk_size=1000000)
        example(reads)

    assert True
