import bionumpy as bnp
import numpy as np
from bionumpy.delimited_buffers import get_bufferclass_for_datatype
from bionumpy.encodings import ACTGEncoding
from bionumpy.groupby import groupby
import matplotlib.pyplot as plt

# Boolean indexing
intervals = bnp.open("example_data/small_interval.bed").read()
print(intervals)
mask = (intervals.end-intervals.start) > 20
print(intervals[mask])

# indexin-slicing
print(intervals[10:-10:2])

# list indexing
reference_sequence = bnp.open("example_data/small_genome.fa.fai")["1"]
indices = np.flatnonzero(reference_sequence[:-1]=="G")
print(indices)
next_letters = reference_sequence[indices+1]
print(next_letters)

np.bincount(bnp.as_encoded_sequence_array(next_letters, ACTGEncoding))

# broadcasting
one_hot_encoding = reference_sequence[:, np.newaxis] == "ACTG"
print(one_hot_encoding)


reads = bnp.open("example_data/big.fq.gz").read()
print(reads)

# reductions
g_count = (reads.sequence == "G").mean(axis=-1)
plt.hist(g_count); plt.show()

# searchsorted
variants = bnp.open("example_data/few_variants.vcf").read()
variants = dict(groupby(variants, "chromosome"))
intervals = dict(groupby(intervals, "chromosome"))
print(variants["1"])

# print(np.where(reference_sequence == "T", bnp.as_sequence_array("U"), reference_sequence))
# where

# broadcasting
# reductions
