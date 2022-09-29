"""
This scripts performs some of the same actions that the commonly
used tools FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
"""

import numpy as np

import bionumpy as bnp
from bionumpy.npdataclassstream import streamable
import matplotlib.pyplot as plt


# @streamable
def get_base_quality_histogram(reads):
    return bnp.bincount(reads.quality.ravel(), minlength=60)


def plot_base_qualities(reads):
    qualities = get_base_quality_histogram(reads)
    plt.plot(qualities)
    plt.show()


@streamable()
def get_gc_content(reads):
    sequences = reads.sequence
    mask = (sequences == "G") | (sequences == "C")
    return np.mean(mask, axis=-1)


def plot_gc_content(reads):
    histogram, _ = bnp.histogram(get_gc_content(reads), bins=50, range=(0, 1))
    plt.plot(histogram)
    plt.show()


@streamable()
def get_quality_scores_as_matrix(reads, limit_at_n_bases=150):
    return reads.quality.as_padded_matrix(side="right", fill_value=0)[:,0:limit_at_n_bases]


def plot_averege_quality_scores_per_base(reads):
    scores = bnp.mean(get_quality_scores_as_matrix(reads), axis=0)
    plt.plot(scores)
    plt.show()


def test():
    examples = [plot_base_qualities, plot_gc_content, plot_averege_quality_scores_per_base]
    for example in examples:
        reads = bnp.open("example_data/big.fq.gz").read_chunks(chunk_size=1000000)
        example(reads)

    assert True


test()
