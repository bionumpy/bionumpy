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

Sequences
~~~~~~~~~
`bnp.Sequences` objects are numpy arrays structured in a way that allows them to hold many sequences of unequal lenght. Under the hood they are `npstructures.RaggedArray` objects that hold byte-arrays, with an encoding that specifies which characters each byte represents. Most of the time, it is not necessary to think about these inner workings, but one can think of them as lists of strings (with the possibility of performing numpy functions on them). The most common way to get `Sequences` objects is to read a file, but they can also be created from lists of strings using the `bnp.as_sequence_array` function::

    >>> bnp.as_sequence_array(["acgttgta", "gcttca", "gttattc"])
    Sequences(acgttgta, gcttca, gttattc)

Encodings
~~~~~~~~~
The main point of bioNumpy is to leverage the computational power of numpy for biological data. A key element to this is to represent different types of biological data as numbers (in numpy arrays). The basic mechanism for doing this is by Encoding classes that can encode data types into numbers, and decode them back in to the data type. A driving idea in bionumpy is to make this happen automatically under the hood, so that a user can choose to ignore the inner workings and (in most cases) relate to sequence data sets in the same way as one would with standard numerical numpy data structures.


Summarization
~~~~~~~~~~~~~
A key application is to extract features from sequence data sets. A large set of interesting features can be computed as functions from sequences to scalar values. Examples are kmer-hashing (kmer->hash-value), minimizers(window->hash-value), string/motif-matching (sequence->bool), Position Weight Matrix scores (sequence->float). bioNumpy provides functionality to apply such functions to rolling windows across large sequence sets, through the `RollableFunction` class. By specifying a broadcastable function in the `__call__` method, the `rolling_window` method will apply the function to all windows in a sequence set. Take the `PositionWeightMatrix` class for instance::


    class PositionWeightMatrix(RollableFunction):
        def __init__(self, matrix, encoding=ACTGEncoding):
            self._matrix = np.asanyarray(matrix)
            self.window_size = self._matrix.shape[-1]
            self._indices = np.arange(self.window_size)
            self._encoding = ACTGEncoding
    
        def __call__(self, sequence: Sequence) -> float:
            sequence = as_sequence_array(sequence, self._encoding)
            scores = self._matrix[sequence, self._indices]
            return scores.sum(axis=-1)

It's `__call__` method specifies how to calculate the score of a sequence. Calling it's rolling_window function will calculate the scores for all windows in a data set::

    >>> import numpy as np
    >>> sequences = bnp.as_sequence_array(["acgttgta", "gcttca", "gttattc"], encoding=bnp.encodings.ACTGEncoding)
    >>> 
    >>> matrix = np.log([[0.1, 0.2],
    ...                  [0.2, 0.3],
    ...                  [0.4, 0.1],
    ...                  [0.3, 0.4]])
    >>> pwm = bnp.position_weight_matrix.PositionWeightMatrix(matrix)
    >>> pwm("ac")
    -3.506557897319982
    >>> pwm(["ac", "cg"])
    array([-3.5065579 , -2.52572864])
    >>> pwm.rolling_window(sequences)
    RaggedArray([[-3.506557897319982, -2.525728644308255, -3.506557897319982, -3.2188758248682006, -1.83258146374831, -3.506557897319982, -2.525728644308255], [-2.4079456086518722, -3.9120230054281455, -3.2188758248682006, -2.120263536200091, -3.2188758248682006], [-3.506557897319982, -3.2188758248682006, -2.525728644308255, -4.605170185988091, -3.2188758248682006, -2.120263536200091]])

Further processing can be done with numpy functions, for instance finding the max score for each sequence in the set::

    >>> pwm.rolling_window(sequences).max(axis=-1)
    array([-1.83258146, -2.12026354, -2.12026354])


kmers
~~~~~
Another example of this concept is the kmer hashing class::

    class KmerEncoding(RollableFunction):
    
        def __init__(self, k, alphabet_size=4):
            self.window_size = k
            self._k = k
            self._alphabet_size = alphabet_size
            self._convolution = self._alphabet_size ** np.arange(self._k)
    
        def __call__(self, sequence: Sequence) -> np.ndarray:
            return sequence.dot(self._convolution)

Here, the `__call__` function specifies how to hash a kmer into a single number. Calling its `rolling_window` method will hash all the kmers in a sequence set.

    >>> bnp.KmerEncoding(3).rolling_window(sequences)
    RaggedArray([[52, 45, 43, 58, 46, 11], [39, 41, 26, 6], [43, 10, 34, 40, 26]])

To count all the 3-mers in the 'reads.fq' sequences we can do as follows:

    >>> fastq_entries_stream = bnp.open("example_data/reads.fq")
    >>> counts = np.zeros(4**3, dtype=int)
    >>> kmer_encoding = bnp.KmerEncoding(3)
    >>> for fastq_entries in fastq_entries_stream:
    ...     kmer_hashes = kmer_encoding.rolling_window(fastq_entries.sequence)
    ...     counts += np.bincount(kmer_hashes.ravel(), minlength=4**3)
    ... 
    >>> counts
    array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0])



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
