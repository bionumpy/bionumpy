Simple quality checking of fastq files
---------------------------------------


In this tutorial we will perform simple quality checking of reads from a fastq file, similarly to what the popular tool FastQC. In addition to BioNumPy, you will also need matplotlib to do some plotting.

We assume you have already followed the introduction part of reading files (see :ref:`reading_files`).

We start by importing all we need:

    >>> import numpy as np
    >>> import bionumpy as bnp
    >>> from bionumpy.npdataclassstream import streamable
    >>> import matplotlib.pyplot as plt


We will be using the `big.fq.gz` file in the example_data folder, but feel free to use any fastq file you like.

The first step is to read our data as chunks:

    >>> reads = bnp.open("example_data/big.fq.gz").read_chunks(chunk_size=1000000)

Note that we now only have a generator object that will give us chunks when we start iterating over it. No data has been read yet.


===========
GC-content
===========

We want to get the GC-content (as a ratio between 0 and 1) within each sequence in a fastq file, and plot a histogram of these numbers.

For each chunk we read from the file, we get the sequences as a RaggedArray where each row is a sequence. Creating a boolean mask of where we have G or C is then as simple as:

    >>> mask = (chunk.sequence == "G") | (chunk.sequence == "C")

`mask` is now still a ragged array with 1 where we have a C or G and 0 elsewhere.

Getting the GC-content for each read can now be done by taking the mean across the reads (last axis) of this mask:

    >>> gc_contents = np.mean(mask, axis=-1)

If we want to do this across all all sequence chunks, we can create a function that does what we want on one chunk and add the streamable decorator:


    >>> @streamable()
    >>> def get_gc_content(reads):
    >>>     sequences = reads.sequence
    >>>     mask = (sequences == ord("G")) | (sequences == ord("C"))
    >>>     return np.mean(mask, axis=-1)

We want to create a histogram of the gc-content values from all chunks. We could call get_gc_content on each chunk, add the results to a list and create a historgram from the final list, but BioNumPy also provides a utility function for creating a histogram from the results from multiple chunks:

    >>> histogram, _ = bnp.histogram(get_gc_content(reads), bins=50, range=(0, 1))
    >>> plt.plot(histogram)
    >>> plt.show()

There is some "magic" happening here that might be useful to understand:

* the @streamable decorator lets us call the get_gc_content on multiple chunks (the result from read_chunks). All it does is to provide an iterable over the results from calling get_gc_content on each chunk.
* the `bnp.histogram` function can take such an iterable and combines the results
* After this code has run, we have iterated over all the chunks in the file, and we need to open the file again and read chunks again if we want to anything else with the file


============================
Histogram of base qualities
============================
If we want to plot a histogram of all the base qualities in all reads, we can use the builtin `bnp.bincount` function. This function does a numpy bincount on each chunk and combines the results.

    >>> base_quality_bincount = bnp.bincount(reads.quality.ravel(), minlength=60)
    >>> plt.plot(base_quality_bincount)
    >>> plt.show()

==============================
Average base quality per base
==============================
In the GC content histogram example, we saw that we can take the mean the rows (axis=-1). If we instead want to find the average base quality for each position in the reads, we can take the mean across the columns (axis=0). Since the reads may have different lengths, we create a padded matrix filled with zeroes. Note that this means that the average base quality is "wrong" after the minimum read length.

    >>> @streamable()
    >>> def get_quality_scores_as_matrix(reads, limit_at_n_bases=150):
    >>>     return reads.quality.as_padded_matrix(side="right", fill_value=0)[:,0:limit_at_n_bases]


    >>> scores = bnp.mean(get_quality_scores_as_matrix(reads), axis=0)
    >>> plt.plot(scores)
    >>> plt.show()

Remember to change the limit_at_n_bases depending on your minimum read length (or how much of the reads you want to plot).

