import numpy as np
from bionumpy.genomic_data.genome import Genome
import bionumpy as bnp


def test(stream=True):
    # read a genome
    genome = Genome.from_file("example_data/hg38.chrom.sizes")

    # read peaks as intervals
    peaks = genome.read_intervals("example_data/ctcf_chr21-22.bed.gz", stream=stream)

    # Get reads as intervals
    reads = genome.read_intervals("example_data/ctcf_chr21-22.bam", stream=stream)

    # reads are just intervals. We can get a pileup across the genome:
    read_pileup = reads.get_pileup()
    print(read_pileup)

    # This pileup can be index, e.g. by our peaks to get a pileup for each peak
    peaks_pileup = read_pileup[peaks]

    # This Pileup can be indexed like a RaggedArray, and we can e.g. get the max pileup
    # value for each peak
    max_pileup_value_per_peak = np.max(peaks_pileup, axis=1)

    # We can filter the peak intervals
    high_peaks = peaks[max_pileup_value_per_peak > 4]
    print(high_peaks)

