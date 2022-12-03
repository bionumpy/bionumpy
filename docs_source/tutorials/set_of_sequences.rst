Simulating sequence datasets
--------

Simulating sequences is very straightforward in bionumpy. Since a sequence arrays have underlying numeric representations that are easy to relate to, one can directly simulate such an underlying representation using any numeric simulation procedure. In addition, bionumpy provides convenience functions that allows to directly simulate sequence data without having to think about any underlying numeric representation. This is the focus of the current tutorial. The function simulate_sequences is an easy way to simulate a set of sequences, by simply specifying a desired alphabet and a dictionary with desired sequence ids as keys and the desired length of each such sequence as value (here simulating 20 sequences with length 10..30):

    >>> import numpy as np
    >>> from bionumpy.simulate import simulate_sequences
    >>> from bionumpy.sequence import match_string
    >>> named_seqs = simulate_sequences('ACGT', {f's{i}':10+i for i in range(20)})
    >>> named_seqs

One can now easily do a variety of analyses on these simulated sequences, e.g. compute the GC content per simulated sequence:
    >>> seqs = named_seqs.sequence
    >>> gc_content_per_seq = np.mean((seqs=='C')|(seqs=='G'), axis=1)
    >>> gc_content_per_seq

If desired, such computed values per sequence can easily be added back as an additional column of the bionumpy data structure:
    >>> named_seqs = named_seqs.add_fields({'gc':gc_content_per_seq}, {'gc':float})
    >>> named_seqs

We can also easily apply a variety of built-in bionumpy functionality on our simulated sequences:
    >>> ac_hits = match_string("AC", seqs)
    >>> ac_hit_sums = np.sum(ac_hits,axis=1)
    >>> ac_hit_sums
    >>> named_seqs = named_seqs.add_fields({'ac_hits':ac_hit_sums}, {'ac_hits':int})
    >>> named_seqs