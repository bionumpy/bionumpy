BioNumPy
======================================

BioNumPy is a Python library for easy and efficient representation and analysis of biological data.
Since BioNumPy builds on the interface of NumPy, people already used to NumPy or array programming should find BioNumPy very easy to get started with.

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

    .. grid-item-card:: :material-regular:`library_books;3em`
        :text-align: center
        :link: topics/sequence_analysis.html

        **Sequence analysis**

        Reading and analysing DNA and protein sequences

    .. grid-item-card::  :material-regular:`location_searching;3em`
        :text-align: center
        :link: topics/kmers.html

        **Kmers**

        Analysing sequence patterns such as kmers, minimzers and motifs

    .. grid-item-card::  :material-regular:`calculate;3em`
        :text-align: center
        :link: topics/genome_arithmetics.html

        **Genome arithmetics**

        Analysing genomic tracks (BED-files, VCFs, GFFs, etc)


.. grid:: 3

    .. grid-item-card:: :material-regular:`hub;3em`
        :text-align: center
        :link: topics/multiomics.html

        **Multiomics**

        Combining data-sets from multiple sources/domains

    .. grid-item-card:: :material-regular:`rocket_launch;3em`
        :text-align: center
        :link: topics/gpu.html

        **GPU-acceleration**

        Ultra-fast sequence analysis using GPU


    .. grid-item-card::  :material-regular:`construction;3em`
        :text-align: center
        :link: topics/extending_bionumpy.html

        **Build on BioNumpy**

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
   source/encoding.rst
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
   developer_guide/encodings.rst
