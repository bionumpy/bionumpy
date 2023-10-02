"""
This scripts performs some of the same actions that the commonly
used tools FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
"""

import numpy as np
import bionumpy as bnp


@bnp.streamable()
def get_base_qualities(reads):
    return reads.quality.ravel()


def plot_base_qualities(reads):
    qualities = bnp.bincount(get_base_qualities(reads), minlength=60)
    return qualities


@bnp.streamable()
def get_gc_content(reads):
    sequences = reads.sequence
    mask = (sequences == "G") | (sequences == "C")
    return np.mean(mask, axis=-1)


def plot_gc_content(reads):
    histogram, _ = bnp.histogram(get_gc_content(reads), bins=50, range=(0, 1))
    return histogram


@bnp.streamable()
def get_quality_scores(reads):
    return reads.quality


def plot_averege_quality_scores_per_base(reads):
    scores = bnp.mean(get_quality_scores(reads), axis=0)
    return scores


def test(do_plot=False):
    examples = [plot_base_qualities, plot_gc_content, plot_averege_quality_scores_per_base]
    for example in examples:
        reads = bnp.open("example_data/big.fq.gz").read_chunks(min_chunk_size=1000000)
        result = example(reads)
        if do_plot:
            import matplotlib.pyplot as plt
            plt.plot(result)
            plt.show()

    assert True


if __name__ == "__main__":
    test(do_plot=True)
