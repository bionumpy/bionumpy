.. _fastqc_tutorial

FastQC-like quality-checking of FASTQ files
------------------------------------------------


In this tutorial we will perform simple quality checking of reads from a fastq file, similarly to what the popular tool FastQC. In addition to BioNumPy, you will also need matplotlib to do some plotting.

We assume you have already followed the introduction part of reading files (see :ref:`reading_files`).

We start by importing all we need:

    >>> import numpy as np
    >>> import bionumpy as bnp
    >>> import matplotlib.pyplot as plt # doctest: +SKIP


We will be using the `big.fq.gz` file in the example_data folder, but feel free to use any fastq file you like.

The first step is to read our data as chunks:

    >>> reads = bnp.open("example_data/big.fq.gz").read_chunks()

Note that we now only have a generator object that will give us chunks when we start iterating over it. No data has been read yet.


===========
GC-content
===========

We want to get the GC-content (as a ratio between 0 and 1) within each sequence in a fastq file, and plot a histogram of these numbers.

For each chunk we read from the file, we get the sequences as a RaggedArray where each row is a sequence. Creating a boolean mask of where we have G or C is then as simple as:

    >>> mask = (chunk.sequence == "G") | (chunk.sequence == "C") # doctest: +SKIP

`mask` is now still a ragged array with 1 where we have a C or G and 0 elsewhere.

Getting the GC-content for each read can now be done by taking the mean across the reads (last axis) of this mask:

    >>> gc_contents = np.mean(mask, axis=-1) # doctest: +SKIP

If we want to do this across all all sequence chunks, we can create a function that does what we want on one chunk and add the streamable decorator:


    >>> @bnp.streamable()
    ... def get_gc_content(reads):
    ...     sequences = reads.sequence
    ...     mask = (sequences == "G") | (sequences == "C")
    ...     return np.mean(mask, axis=-1)

We want to create a histogram of the gc-content values from all chunks. We could call get_gc_content on each chunk, add the results to a list and create a histogram from the final list, but BioNumPy also provides a utility function for creating a histogram from the results from multiple chunks:

    >>> histogram, _ = bnp.histogram(get_gc_content(reads), bins=50, range=(0, 1))
    >>> plt.plot(histogram) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

There is some "magic" happening here that might be useful to understand:

* the @streamable decorator lets us call the get_gc_content on multiple chunks (the result from read_chunks). All it does is to provide an iterable over the results from calling get_gc_content on each chunk.
* the `bnp.histogram` function can take such an iterable and combines the results
* After this code has run, we have iterated over all the chunks in the file, and we need to open the file again and read chunks again if we want to anything else with the file


============================
Histogram of base qualities
============================
If we want to plot a histogram of all the base qualities in all reads, we can use the builtin `bnp.bincount` function. This function does a numpy bincount on each chunk and combines the results.

    >>> @bnp.streamable()
    ... def get_base_qualities(reads):
    ...     return reads.quality.ravel()
    >>> reads = bnp.open("example_data/big.fq.gz").read_chunks()
    >>> base_quality_bincount = bnp.bincount(get_base_qualities(reads), minlength=60)
    >>> plt.plot(base_quality_bincount) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

==============================
Average base quality per base
==============================
In the GC content histogram example, we saw that we can take the mean the rows (axis=-1). If we instead want to find the average base quality for each position in the reads, we can take the mean across the columns (axis=0).

.. code-block:: python

    @bnp.streamable()
    def get_quality_scores(reads):
        return reads.quality

    reads = bnp.open("example_data/big.fq.gz").read_chunks()
    scores = bnp.mean(get_quality_scores(reads), axis=0)
    print(scores)
    plt.plot(scores)
    plt.show()


Remember to change the limit_at_n_bases depending on your minimum read length (or how much of the reads you want to plot).

