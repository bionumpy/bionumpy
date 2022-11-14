.. _topic_kmers


Using BioNumPy to analyse kmers, minimzers and motifs
========================================================

BioNumPy has functionality for easily finding kmers and minimizers in a set of sequences. For instance, finding all 3-mers and counting the number of `ACT` in some sequences is as simple as:

    >>> import bionumpy as bnp
    >>> sequences = bnp.as_encoded_array(["ACTG", "GGGACT", "G"], bnp.DNAEncoding)
    >>> kmers = bnp.kmers.get_kmers(sequences, 3)
    >>> kmers
    encoded_ragged_array([[ACT, CTG],
                          [GGG, GGA, GAC, ACT],
                          []], 3merEncoding(AlphabetEncoding('ACGT')))
    >>> counts = bnp.count_encoded(kmers, axis=None)
    >>> counts["ACT"]
    array(2)

The same functionality works for minimizers:

    >>> bnp.minimizers.get_minimizers(sequences, k=2, window_size=4)
    encoded_ragged_array([[AC],
                          [GA, GA, GA],
                          []], 2merEncoding(AlphabetEncoding('ACGT')))



Analysing kmers and minimizers in big data sets
-------------------------------------------------
BioNumPy is able to quite efficiently find all kmers or minimizers in large sequence data sets. Since the resulting kmers are stored numerically in NumPy arrays, NumPy functionality can be used to analyse the resulting kmers/minimizers.

Below in an example of finding all 31-mers in a fastq file and ....




[What to read next]

