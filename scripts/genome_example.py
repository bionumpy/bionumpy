import numpy as np
from bionumpy.genomic_data.genome import Genome
import bionumpy as bnp


def test():
    # read a genome
    genome = Genome.from_file("example_data/hg38.chrom.sizes")
    print(genome)

    # read a bed-file as a mask with this genome
    peak_mask = genome.read_intervals("example_data/ctcf.bed.gz").get_mask()
    print(peak_mask)

    # We can get the interval dataclass for this GenomicTrack
    data = peak_mask.get_data()
    print("Average inteval length: ", np.mean(data.stop-data.start))

    # GenomicTracks can be indexed
    chr1_mask = peak_mask["chr1"]
    print(chr1_mask[0:10])

    # GenomicTracks can be indexes with other GenomicTracks
    ctcf = genome.read_intervals("example_data/ctcf.bed.gz").get_mask()
    znf = genome.read_intervals("example_data/znf263.bed.gz")

    print(type(ctcf[znf]))
    subset = ctcf[znf]

