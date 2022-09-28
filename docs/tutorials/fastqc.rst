Simple quality checking of fastq files
---------------------------------------


In this tutorial we will perform simple quality checking of reads from a fastq file, similarly to what the popular tool FastQC. In addition to BioNumPy, you will also need matplotlib to do some plotting.

We start by importing all we need:

    import numpy as np
    import bionumpy as bnp
    from bionumpy.npdataclassstream import streamable
    import matplotlib.pyplot as plt


We will be using the big.fq.gz file in the example_data folder, but feel free to use any fastq file you like.

The first stepis to read our data as chunks:

    reads = bnp.open("example_data/big.fq.gz").read_chunks(chunk_size=1000000)

Note that we now only have a generator object that will give us chunks when we start iterating over it. No data has been read yet.

You can iterate over the chunks using a simple for-loop:



GC-content
-----------
Computing the GC-content from a single sequences is just to count the number of Cs and Gs and dividing by the sequence length.

For each chunk in BioNumPy, we have sequences represented as a RaggedArray which can be thought of a matrix where each row may have different length. Each row is in our case a sequence. This datastructure supports many of the same operations that a normal NumPy matrix supports, i.e. masking out elements:

```
for chunk in reads:
    sequences = chunk.sequence
    print(sequences == ord("C"))
```

If you try to run the above code,



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
    mask = (sequences == ord("G")) | (sequences == ord("C"))
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
