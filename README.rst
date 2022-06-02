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
The main function for file input and output is the `bnp_open` function, which by default deduces the file type based on the suffix, and creates a buffered reader of a sensible type. For instance::

    >>> 

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
