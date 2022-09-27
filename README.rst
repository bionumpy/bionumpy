========
BioNumPy
========


.. image:: https://img.shields.io/pypi/v/bionumpy.svg
        :target: https://pypi.python.org/pypi/bionumpy

.. image:: https://github.com/knutdrand/bionumpy/actions/workflows/python-install-and-test.yml/badge.svg
        :target: https://github.com/knutdrand/bionumpy/actions/
        :alt: Build and test status

.. image:: https://readthedocs.org/projects/bionumpy/badge/?version=latest
        :target: https://bionumpy.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status


What is BioNumPy?
----------------
BioNumPy is a toolkit, built on top of NumPy, for enabling array programming on biological data in Python. BioNumPy aims to make it easy to read common bioinformatics file formats efficiently into NumPy-like data structures that enable efficient operations and analysis of the data. Working in BioNumPy should feel much like working in NumPy.

.. code-block:: python
    a = "hei"
    print(a)


Why BioNumPy?
-------------
* There are no existing packages in Python for getting biological data sets efficiently into NumPy-like data structures
* Current packages for working with biological data do not use NumPy in an efficient way (e.g. individual sequences are stored as separate NumPy arrays, not together in shared arrays).


Getting started
----------------

1. Install with pip:

>>> pip install bionumpy

2. Read the introduction below
3. Check out the more advanced examples in the documentation


Overview
------------
Here is a brief overview of what you need to know in order to get started. We refer to the documentation for more details.

The features of BioNumPy can roughly be divided into two:

1. Reading biological data sets (e.g. fasta, vcf, bed) into NumPy-like objects
2. Analysing, modifying, filtering these NumPy-like objects in meaningful ways

BioNumPy also supports writing most data types to file.


Reading and writing files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The main function for file input and output is the `bnp.open` function, which by default deduces the file type based on the suffix, and creates a buffered reader of a sensible type. For instance::

    >>> import bionumpy as bnp
    >>> 
    >>> file = bnp.open("example_data/reads.fq")
    >>> for fastq_entries in file.get_chunks():
    ...     print(fastq_entries)
    ... 
    SequenceEntryWithQuality(name=Sequences(headerishere, anotherheader), sequence=Sequences(CTTGTTGA, CGG), quality=Sequences(!!!!!!!!, ~~~))

A couple of things to note:

* the `bnp.open` function understands that the `.fq` suffix means it's a fastq file and reads it as such.
* the method `get_chunk()` on the File object will return chunks from the file (in this particular case, the file is so small that we only get one chunk).
* the data from each chunk is delivered as an `npdataclass` where there is one array-like object for each field. This means that if we only care about the acutal DNA-sequences, we can operate on that on it's own.
* in the above example, we could get the sequences as an NumPy-like object (a RaggedArray) like this: `fastq_entries.sequence`)

Another example would be to read in a vcf-file::

    >>> file = bnp.open("example_data/variants.vcf")
    >>> for chromosome, variants in file.get_chunks():
    ...     print(variants)
    ... 
    Variant(chromosome=['chr1' 'chr1' 'chr1'], position=array([883624, 887559, 887800]), ref_seq=Sequences(A, A, A), alt_seq=Sequences(G, C, G))
    Variant(chromosome=['chr9' 'chr9' 'chr9'], position=array([883624, 887559, 887800]), ref_seq=Sequences(A, A, A), alt_seq=Sequences(G, C, G))


Again, `bnp.open` recognizes that this is a vcf file, and again it chooses an appropriate format to output it in.

Sequences
~~~~~~~~~
`bnp.Sequences` objects are numpy arrays structured in a way that allows them to hold many sequences of unequal lenght. Under the hood they are `npstructures.RaggedArray` objects that hold byte-arrays, with an encoding that specifies which characters each byte represents. Most of the time, it is not necessary to think about these inner workings, but one can think of them as lists of strings (with the possibility of performing numpy functions on them). The most common way to get `Sequences` objects is to read a file, but they can also be created from lists of strings using the `bnp.as_sequence_array` function::

    >>> bnp.as_sequence_array(["acgttgta", "gcttca", "gttattc"])
    Sequences(acgttgta, gcttca, gttattc)



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
