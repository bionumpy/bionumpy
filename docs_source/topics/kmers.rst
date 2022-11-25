.. _topic_kmers


Using BioNumPy to analyse kmers, minimzers and motifs
========================================================

BioNumPy has functionality for easily finding kmers and minimizers in a set of sequences. For instance, finding all 3-mers and counting the number of `ACT` in some sequences is as simple as:

import bionumpy.encoded_array_functions    >>> import bionumpy as bnp
    >>> sequences = bnp.as_encoded_array(["ACTG", "GGGACT", "G"], bnp.DNAEncoding)
    >>> kmers = bnp.sequence.get_kmers(sequences, 3)
    >>> kmers
    encoded_ragged_array([[ACT, CTG],
                          [GGG, GGA, GAC, ACT],
                          []], 3merEncoding(AlphabetEncoding('ACGT')))
    >>> counts = bnp.count_encoded(kmers, axis=None)
    >>> counts["ACT"]
    array(2)

The same functionality works for minimizers:

    >>> bnp.sequence.get_minimizers(sequences, k=2, window_size=4)
    encoded_ragged_array([[AC],
                          [GA, GA, GA],
                          []], 2merEncoding(AlphabetEncoding('ACGT')))

Note that the kmer hashes and minimizers are numerical internally, and how they are computed will depend on the encoding of the sequences. If your sequences is encoded with an encoding with alphabet size **4** (e.g.

import bionumpy.encoded_array_functions    >>> import bionumpy as bnp
    >>> sequences = bnp.as_encoded_array(["ACTG", "GGGACT", "G"], bnp.DNAEncoding)
    >>> kmers = bnp.sequence.get_kmers(sequences, 3)
    >>> kmers
    encoded_ragged_array([[ACT, CTG],
                          [GGG, GGA, GAC, ACT],
                          []], 3merEncoding(AlphabetEncoding('ACGT')))
    >>> counts = bnp.count_encoded(kmers, axis=None)
    >>> counts["ACT"]
    array(2)

The same functionality works for minimizers:

    >>> bnp.sequence.get_minimizers(sequences, k=2, window_size=4)
    encoded_ragged_array([[AC],
                          [GA, GA, GA],
                          []], 2merEncoding(AlphabetEncoding('ACGT')))

Note that the kmer hashes and minimizers are numerical internally, and how they are computed will depend on the encoding of the sequences. If your sequences is encoded with an encoding with alphabet size **4** (e.g.

import bionumpy.encoded_array_functions    >>> import bionumpy as bnp
    >>> sequences = bnp.as_encoded_array(["ACTG", "GGGACT", "G"], bnp.DNAEncoding)
    >>> kmers = bnp.sequence.get_kmers(sequences, 3)
    >>> kmers
    encoded_ragged_array([[ACT, CTG],
                          [GGG, GGA, GAC, ACT],
                          []], 3merEncoding(AlphabetEncoding('ACGT')))
    >>> counts = bnp.count_encoded(kmers, axis=None)
    >>> counts["ACT"]
    array(2)

The same functionality works for minimizers:

    >>> bnp.sequence.get_minimizers(sequences, k=2, window_size=4)
    encoded_ragged_array([[AC],
                          [GA, GA, GA],
                          []], 2merEncoding(AlphabetEncoding('ACGT')))

Note that the kmer hashes and minimizers are numerical internally, and how they are computed will depend on the encoding of the sequences. If your sequences is encoded with an encoding with alphabet size **4** (e.g.

import bionumpy.encoded_array_functions    >>> import bionumpy as bnp
    >>> sequences = bnp.as_encoded_array(["ACTG", "GGGACT", "G"], bnp.DNAEncoding)
    >>> kmers = bnp.sequence.get_kmers(sequences, 3)
    >>> kmers
    encoded_ragged_array([[ACT, CTG],
                          [GGG, GGA, GAC, ACT],
                          []], 3merEncoding(AlphabetEncoding('ACGT')))
    >>> counts = bnp.count_encoded(kmers, axis=None)
    >>> counts["ACT"]
    array(2)

The same functionality works for minimizers:

    >>> bnp.sequence.get_minimizers(sequences, k=2, window_size=4)
    encoded_ragged_array([[AC],
                          [GA, GA, GA],
                          []], 2merEncoding(AlphabetEncoding('ACGT')))

Note that the kmer hashes and minimizers are numerical internally, and how they are computed will depend on the encoding of the sequences. If your sequences is encoded with an encoding with alphabet size **4** (e.g.

import bionumpy.encoded_array_functions    >>> import bionumpy as bnp
    >>> sequences = bionumpy.encoded_array_functions.as_encoded_array(["ACTG", "GGGACT", "G"], bnp.DNAEncoding)
    >>> kmers = bnp.sequence.get_kmers(sequences, 3)
    >>> kmers
    encoded_ragged_array([[ACT, CTG],
                          [GGG, GGA, GAC, ACT],
                          []], 3merEncoding(AlphabetEncoding('ACGT')))
    >>> counts = bnp.count_encoded(kmers, axis=None)
    >>> counts["ACT"]
    array(2)

The same functionality works for minimizers:

    >>> bnp.sequence.get_minimizers(sequences, k=2, window_size=4)
    encoded_ragged_array([[AC],
                          [GA, GA, GA],
                          []], 2merEncoding(AlphabetEncoding('ACGT')))

Note that the kmer hashes and minimizers are numerical internally, and how they are computed will depend on the encoding of the sequences. If your sequences is encoded with an encoding with alphabet size **4** (e.g.

import bionumpy.encoded_array_functions    >>> import bionumpy as bnp
    >>> sequences = bionumpy.encoded_array_functions.as_encoded_array(["ACTG", "GGGACT", "G"], bnp.DNAEncoding)
    >>> kmers = bnp.sequence.get_kmers(sequences, 3)
    >>> kmers
    encoded_ragged_array([[ACT, CTG],
                          [GGG, GGA, GAC, ACT],
                          []], 3merEncoding(AlphabetEncoding('ACGT')))
    >>> counts = bnp.count_encoded(kmers, axis=None)
    >>> counts["ACT"]
    array(2)

The same functionality works for minimizers:

    >>> bnp.sequence.get_minimizers(sequences, k=2, window_size=4)
    encoded_ragged_array([[AC],
                          [GA, GA, GA],
                          []], 2merEncoding(AlphabetEncoding('ACGT')))

Note that the kmer hashes and minimizers are numerical internally, and how they are computed will depend on the encoding of the sequences. If your sequences is encoded with an encoding with alphabet size **4** (e.g.

    >>> import bionumpy as bnp
    >>> sequences = bnp.as_encoded_array(["ACTG", "GGGACT", "G"], bnp.DNAEncoding)
    >>> kmers = bnp.sequence.get_kmers(sequences, 3)
    >>> kmers
    encoded_ragged_array([[ACT, CTG],
                          [GGG, GGA, GAC, ACT],
                          []], 3merEncoding(AlphabetEncoding('ACGT')))
    >>> counts = bnp.count_encoded(kmers, axis=None)
    >>> counts["ACT"]
    array(2)

The same functionality works for minimizers:

    >>> bnp.sequence.get_minimizers(sequences, k=2, window_size=4)
    encoded_ragged_array([[AC],
                          [GA, GA, GA],
                          []], 2merEncoding(AlphabetEncoding('ACGT')))

Note that the kmer hashes and minimizers are numerical internally, and how they are computed will depend on the encoding of the sequences. If your sequences is encoded with an encoding with alphabet size **4** (e.g. `bnp.DNAEncoding`) `get_kmers` is optimized to run a lot faster by performing bit-operations.


Analysing kmers and minimizers in big data sets
-------------------------------------------------
BioNumPy is able to quite efficiently find all kmers or minimizers in large sequence data sets. Since the resulting kmers are stored numerically in NumPy arrays, NumPy functionality can be used to do analysis.

Below in an example of finding all 31-mers in a fastq file and



Read next
----------

[What to read next]

