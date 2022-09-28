.. _reading_files:

Reading files
---------------

There are three main ways of getting your data into memory in BioNumPy. Which way you choose mainly depends on how big your data is:

Method 1: Reading all your data at once
========================================
If you are working with data sets that are small enough to fit into memory, you can read the whole file. The benefit of this approach is that you don't need to take into account that the data has been split into chunks (next methods):

    >>> data = open("filename").read()
    >>> print(data)

In the above example, `data` will be an NpDataClass object with various fields, depending on the file type. For instance, if you read a fastq file, you will be able to access the sequences (`data.sequence`) and the bsae qualities (`data.quality`).


Method 2: Reading a single chunk from a large file
===================================================
If you have a large file, and don't want to read all of it, you can read a chunk. This can be useful if you only want to get to know your data or perform an analysis on a small subset of your data:

    >>> file = open("reads.fastq")
    >>> chunk = file.read_chunk(chunk_size=10000000)

Here `chunk_size` is the number of bytes to read. If the file is small enough, you may get the whole file by setting chunk size to a big number, but if you want to make sure you read the whole file as one single cunk you should choose method 1 instead.


Method 3: Reading whole file as chunks
========================================
This is the preferred way you should use for all large files (typically fasta files, fastq files and bam files). The idea is to read the whole file as chunks that are iterated over. This way, only one single chunk will be kept in memory at the time.

    >>> file = bnp.open("file.fastq")
    >>> for chunk in file.read_chunks():
    >>>    # do stuff with this chunk, e.g. take the man
    >>>    # of the base qualities


Remember that when one chunk has been processed, that chunk is "lost" and the generator continues to the next chunk. A typical implementation will thus be to write a function that does all you want to do one one chunk, run that function on each chunk and summarize the results in some way.

Continue to the guide on :ref:`working_with_big_data` to see the recommended way of working with chunks in BioNumPy.

