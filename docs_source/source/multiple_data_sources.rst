.. _intervals:

========================================
Working with Multiple Files/Data Sources
========================================

A lot of analyses requires working on multiple data sources that are annotated on the same reference sequences. Especially when workin with reference genomes. A common example is analysing soms variants from a vcf file, but only those variants that fall within some intervals specified in a bed file. One way to do this in BioNumPy is to use the `MultiStream` class. This is designed to take in multiple data streams or indexed data and deliver synchronized streams that return data corresponding to the same reference sequence/chromosome in each chunk. It is best shown in an example.

    >>> import bionumpy as bnp
    >>> variants = bnp.open("example_data/few_variants.vcf").read_chunks()
    >>> intervals = bnp.open("example_data/small_interval.bed").read_chunks()
    >>> reference = bnp.open_indexed("example_data/small_genome.fa")

These variants and intervals come in chunks, and stream of chunks are not really satisfied. In addition the reference genome is read as an indexed file, and so is able to give us arbitrary parts of the reference genome at demand. We can synch these three sources up using MultiStream:

    >>> multistream = bnp.MultiStream(reference.get_contig_lengths(),
    ...                               sequence=reference,
    ...                               variants=variants,
    ...                               intervals=intervals)


The attributes we specify for the multistream are now synched streams that give data corresponding to the seqeunces listed int `reference.get_sequence_lengths()` one at a time. This means we can give these streams to any function with the `streamable` decorator. For example:

    >>> from bionumpy.intervals import get_boolean_mask
    >>> @bnp.streamable(sum)
    ... def get_letter_counts_around_variants(reference_sequence, variants, intervals):
    ...     mask = get_boolean_mask(intervals, len(reference_sequence))
    ...     positions = variants.position
    ...     positions = positions[mask[positions]] # Only look at variants inside intervals
    ...     surrounding_sequences = reference_sequence[positions-10:positions+10].ravel()
    ...     return bnp.count_encoded(bnp.as_encoded_array(surrounding_sequences, bnp.DNAEncoding))
    
    
    >>> get_letter_counts_around_variants(multistream.sequence, multistream.variants, multistream.intervals)
    EncodedCounts(alphabet=['A', 'C', 'G', 'T'], counts=array([155, 127, 116, 122]), row_names=None)
