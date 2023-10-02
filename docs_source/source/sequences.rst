.. _sequences:

=========
Sequences
=========

`EncodedArray`
==============

Sequence data in BioNumPy are are represented by `EncodedArray` objects. These are basically `numpy` arrays of integers, that have an encoding that specifies which character each integer represent. This representataion allows us to do fast `numpy` operation on the seqeunces, while still allowing for human readable representation of them. The easiest way to create an `EncodedArray` is to use the `as_encoded_array` function.

    >>> import bionumpy as bnp
    >>> encoded_array = bnp.as_encoded_array("actggtcc")
    >>> encoded_array
    encoded_array('actggtcc')

We see that the encoded array represents the text "actggtcc", but under the hood these are just integers with an ASCII encoding. You seldom have to think about this internal representation, but it is what allow us to write for instance:

    >>> print(encoded_array == "g")
    [False False False  True  True False False False]
    
And get numpy-fast performance for the query. We can also use numpy-like indexing on encoded arrays, so that we can for instance trim the first and last two characters from the sequence:

    >>> encoded_array[2:-2]
    encoded_array('tggt')
    
`EncodedRaggedArray`
====================
When working with multiple seqeunces we usually have to use `EncodedRaggedArray` objects. These are much like `EncodedArray` objects, but instead uses `npstructures.RaggedArray` to store the integers. This allows us to store seqeunces of differing lengths. The easiest way to create an `EncodedRaggedArray` is to use the `as_encoded_array` function on a list of strings:

    >>> encoded_ragged_array = bnp.as_encoded_array(["ctt", "actg", "ag"])
    >>> encoded_ragged_array
    encoded_ragged_array(['ctt',
                          'actg',
                          'ag'])

These objects also behave very much like numpy arrays, in indexing and broadcasting. For instance, to get the 2nd character of the first and third seqeunce:

    >>> encoded_ragged_array[[0, 2], 1]
    encoded_array('tg')

Reading sequences from file
===========================
Usually we get sequences directly from file. BioNumPy supports a range of file formats containing sequence data including fasta, fastq, indexed fasta and bam files.


Reading sequence entries (.fa, .fq, .bam)
-----------------------------------------
To read in a set of sequence entries, we can just use the `bnp.open` method (here we read a fastq file, but this works the same for fasta and bam files):

    >>> entries = bnp.open("example_data/reads.fq").read()
    >>> entries
    SequenceEntryWithQuality with 2 entries
                         name                 sequence                  quality
                 headerishere                 CTTGTTGA        [2 2 2 2 2 2 2 2]
                anotherheader                      CGG               [93 93 93]

We see we get all the entries in the file, with the corresponding fields. The `sequence` field here is an `EncodedRaggedArray` and thus supports numpy-like indexing etc:

    >>> (entries.sequence=="T").sum(axis=-1) # Count the number of T's
    array([4, 0])


Reading indexed files (.fa.fai)
-------------------------------
When reading a reference genome, we often can't read in the whole file (using `.read()`) and it doesn't make sense to read in the chromosome-sequences as entries in chunks (using `.read_chunks()`) but we rather want to read specific parts of the genome. In these cases we can use an index to be able to read specific regions. We then use `bnp.open_indexed` function:

    >>> reference_sequence = bnp.open_indexed("example_data/small_genome.fa")
    >>> reference_sequence["2"][10:20]
    encoded_array('atattagcca')

Functions workin on sequences
=============================

A set of functions working on sequences are gathered in the `bnp.seqeunce` module. 
