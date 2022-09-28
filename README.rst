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
    >>> for fastq_entries in file.read_chunks():
    ...     print(fastq_entries)
    ... 
    SequenceEntryWithQuality(name=Sequences(headerishere, anotherheader), sequence=Sequences(CTTGTTGA, CGG), quality=Sequences(!!!!!!!!, ~~~))

A couple of things to note:

* the `bnp.open` function understands that the `.fq` suffix means it's a fastq file and reads it as such.
* the method `read_chunks()` on the File object will return chunks from the file (in this particular case, the file is so small that we only get one chunk).
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


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
