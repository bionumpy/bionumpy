========
bionumpy
========


.. image:: https://img.shields.io/pypi/v/bionumpy.svg
        :target: https://pypi.python.org/pypi/bionumpy

.. image:: https://img.shields.io/travis/knutdrand/bionumpy.svg
        :target: https://travis-ci.com/knutdrand/bionumpy

.. image:: https://readthedocs.org/projects/bionumpy/badge/?version=latest
        :target: https://bionumpy.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Library for working with biological sequence data as numpy arrays.


* Free software: MIT license
* Documentation: https://bionumpy.readthedocs.io.


Features
--------

bioNumpy is a library for handling large biological sequence date sets in Python. The main features of bioNumpy are input and output from common file formats (fastq, fasta, bed, bam, sam...) (IO); Encoding of common data types into numerical arrays (Encoding) and summarizing/counting features in a data sets (kmers, minimizers, gc content, string matching,,,).

IO
~~
The main function for file input and output is the `open` function, which by default deduces the file type based on the suffix, and creates a buffered reader of a sensible type. For instance::

    >>> import bionumpy as bnp
    >>> 
    >>> fastq_entries_stream = bnp.open("example_data/reads.fq")
    >>> for fastq_entries in fastq_entries_stream:
    ...     print(fastq_entries)
    ... 
    SequenceEntryWithQuality(name=Sequences(headerishere, anotherheader), sequence=Sequences(CTTGTTGA, CGG), quality=Sequences(!!!!!!!!, ~~~))

A couple of things to note: Firstly the `bnp.open` function understands that the `.fq` suffix means it's a fastq file and reads it as such. It also knowns that usually a `fastq` file is usually too big to keep in memory so instead of reading in the whole file, it reads in chunks of data and delivers it as a stream (In this particular case, the file is so small that we only get one chunk). Third, it delivers the data for each chunk as a `npdataclass` where there is one array-like object for each field. This means that if we only care about the acutal DNA-sequences, we can operate on that on it's own.

Another example would be to read in a vcf-file as such::

    >>> variants_stream = bnp.open("example_data/variants.vcf")
    >>> for chromosome, variants in variants_stream:
    ...     print(variants)
    ... 
    Variant(chromosome=['chr1' 'chr1' 'chr1'], position=array([883624, 887559, 887800]), ref_seq=Sequences(A, A, A), alt_seq=Sequences(G, C, G))
    Variant(chromosome=['chr9' 'chr9' 'chr9'], position=array([883624, 887559, 887800]), ref_seq=Sequences(A, A, A), alt_seq=Sequences(G, C, G))


Again, `bnp.open` recognizes that this is a vcf file, and again it chooses an appropriate format to output it in. For vcf file this is a stream variant chunks, where each chunk is the variants for one chromosome. This format together with the `ChromosomeMap` decorator makes it easy to work with genomic data.

Encodings
~~~~~~~~~
The main point of bioNumpy is to leverage the computational power of numpy for biological data. A key element to this is to represent different types of biological data as numbers (in numpy arrays). The basic mechanism for doing this is by Encoding classes that can encode data types into numbers, and decode them back in to the data type.


Summarization
~~~~~~~~~~~~~
A key application is to extract features from data sets. A large set of interresting features can be computed as functions from sequences to scalar values. Examples are kmer-hashing (kmer->hash-value), minimizers(window->hash-value), string/motif-matching (sequence->bool), Position Weigh Matrix scores (sequence->float). bioNumpy provides functionality to apply such functions to rolling windows across large  seuqence sets, through the `rolling_window_function` decorator. The `rolling_window_function` decorator applies the function to all windows in a `ndarray` or `RaggedArray` of sequences. The only requirement is that the function is broadcastable::

    sequences = Sequences(["acgttgta", "gcttca", "gttattc"], encoding=ACGTEncoding)
    
    matrix = np.log([[0.1, 0.2],
                     [0.2, 0.3],
                     [0.4, 0.1],
                     [0.3, 0.4]])
    
    pwm = PositionWeightMatrix(matrix)
    pwm.calculate_score("ac")
    pwm.calculate_score(["ac", "cg"])
    rolling_pwm(sequences, 2, pwm)

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
