BioNumPy
======================================

BioNumPy is a Python library for easy and efficient representation and analysis of biological data.
Since BioNumPy builds on the interface of NumPy, people already used to NumPy or array programming should find BioNumPy very easy to get started.

..
    With BioNumPy, our goal is that everyone should be able to write simple, clean code that scales well to large biological data sets.
**Getting started with BioNumpy takes only a minute:**

1) Install:

.. code-block:: bash

    pip install bionumpy

2) Read your data:

    >>> import numpy as np
    >>> import bionumpy as bnp
    >>> reads = bnp.open("example_data/small.fa").read()
    >>> reads
    SequenceEntry with 3 entries
                         name                 sequence
                        read1  ACACATCACAGCTACGACGA...
                        read2  AACACTTGGGGGGGGGGGGG...
                        read3  AACTGGACTAGCGACGTACT...

3) Analyse it:

    >>> gc_content = np.mean((reads.sequence == "C") | (reads.sequence == "G"))
    >>> gc_content
    0.5526315789473685

BioNumpy can be used to analyse a wide range of data. Follow one of the guides below:

What can you do with BioNumpy?
----------------------------------

.. grid:: 3

    .. grid-item-card:: Sequence analysis
        :link: https://example.com

        Reading and analysing DNA and protein sequences

    .. grid-item-card::  Kmers
        :link: topics/kmers.html

        Analysing sequence patterns such as kmers, minimzers and motifs

    .. grid-item-card::  Genome arithmetics
        :link: https://example.com

        Analysing genomic tracks (BED-files, VCFs, GFFs, etc)


.. grid:: 3

    .. grid-item-card::  Multi-omics
        :link: https://example.com

        Combining data-sets from multiple sources/domains

    .. grid-item-card::  GPU-acceleration
        :link: https://example.com

        Ultra-fast sequence analysis using GPU


    .. grid-item-card::  Build on BioNumpy
        :link: https://example.com

        Combine core functionality to support your use-case





If you are interested in learning more, learn more about BioNumPy (todo: link to source/concepts), check out the tutorials (TODO: LINK), how to work with different file formats, or check out the API documentation.


Documentation
=================

.. toctree::
   :maxdepth: 0
   :glob:
   :titlesonly:
   :caption: Tutorials

   tutorials/*

.. toctree::
   :maxdepth: 2
   :glob:
   :caption: BioNumPy Concepts

   source/reading_files.rst
   source/working_with_big_data.rst
   source/supported_file_formats.rst
   source/sequences.rst
   source/intervals.rst
   source/multiple_data_sources.rst
   source/broadcastable_functions.rst
   source/rollable_functions.rst
   source/summarization.rst


.. toctree::
   :maxdepth: 2
   :caption: API documentation
   :glob:

   modules


.. toctree::
   :maxdepth: 1
   :caption: Developer Guide

   developer_guide/getting_started.rst
   developer_guide/setting_up_development_environment.rst
   developer_guide/testing.rst
   developer_guide/making_examples.rst
   developer_guide/writing_documentation.rst
   developer_guide/design_principles.rst
