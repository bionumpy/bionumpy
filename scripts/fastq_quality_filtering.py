import numpy as np
import bionumpy as bnp
from npstructures import SeqArray, RaggedArray


fastq_stream = bnp.open("example_data/reads.fq")

for fastq_reads in fastq_stream:
    # find sequences with mean quality below 20
    mask = np.mean(fastq_reads.quality, axis=-1) < 20
    print(mask)

    # find sequences with lowest quality base below 20
    mask = np.min(fastq_reads.quality, axis=-1) < 100
    print(mask)

    # Plot average sequence quality per read
    avg_qualities = np.mean(fastq_reads.quality, axis=-1)
    print(avg_qualities.encoding)
    print(repr(avg_qualities))
    print(avg_qualities)

    # Sequence length distribution
    #lengths = fastq_reads.sequence.shape
    #print(lengths)

    # todo: Convert ragged array to padded matrix would be nice??

    #print(max_qualities)
    #print(np.max(fastq_reads.quality))
    #keep = np.where(fastq_reads.quality > 0)

    #print(keep)

    #passed = fastq_reads.replace(np.where(fast))