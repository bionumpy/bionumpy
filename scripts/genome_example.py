import numpy as np
from bionumpy.genomic_data.genome import Genome
import bionumpy as bnp


def test():
    # read a genome
    genome = Genome.from_file("example_data/hg38.chrom.sizes")
    print(genome)

    # read peaks as intervals
    peaks = genome.read_intervals("example_data/ctcf_chr21-22.bed.gz")

    # Get reads as intervals
    reads = genome.read_intervals("example_data/ctcf_chr21-22.bam")

    # reads are just intervals. Get a pileup (GenomeArray) across the genome
    pileup = reads.get_pileup()

    # Fetch the pileup in peak areas
    peaks_pileup = pileup[peaks]

    print(peaks_pileup)

    return
    # Get a GenomeArray boolean mask (True where we have peaks)
    mask = peaks.get_mask()

    return
    # We can get the interval dataclass for this GenomicTrack
    data = peak_mask.get_data()
    print("Average inteval length: ", np.mean(data.stop-data.start))

    # GenomicTracks can be indexed
    chr1_mask = peak_mask["chr1"]
    print(chr1_mask[0:10])

    # GenomicTracks can be indexes with Intervals
    ctcf = genome.read_intervals("example_data/ctcf.bed.gz").get_mask()
    znf = genome.read_intervals("example_data/znf263.bed.gz")

    print(type(ctcf))
    print(type(znf))
    #print(type(ctcf[znf]))
    #print(znf[ctcf])

