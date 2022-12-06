.. _supported_file_formats:

Supported file formats
-----------------------------------

This is a list of  currently supported file formats in BioNumPy. Reading files with any of these extensions with `bnp.open` will make BioNumPy automatically detect the file type and read the data into an appropriate data structure (which will be a dataclass-like object with fields).

* vcf
* bed
* fasta / fa
* bed
* fasta / fa
* gfa (limited support only)
* gff
* gtf
* gff3
* sam / bam

=======
Example
=======
We open a bed file, read one chunk and print a description of that chunk:

.. testcode::

    import bionumpy as bnp
    data = bnp.open("example_data/test.bed")
    chunk = data.read_chunk()
    print(chunk)

The above example should work wih any of the supported file formats.

This shows us that we have a a chunk of 71 intervals, and we get to see the first of these:

.. testoutput::

    Interval with 71 entries
                   chromosome                    start                     stop
                           17                  7512371                  7512447
                           17                  7512377                  7512453
                           17                  7512393                  7512469
                           17                  7512420                  7512496
                           17                  7512422                  7512473
                           17                  7512425                  7512501
                           17                  7512428                  7512504
                           17                  7512474                  7512550
                           17                  7512537                  7512613
                           17                  7512559                  7512635


Implementing a new file format
------------------------------
This guide is not written.
