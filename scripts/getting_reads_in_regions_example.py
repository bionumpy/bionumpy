"""
This example shows how to subset alignments in a BED-file
on regions specified in another BED-file. Using the raggedslice-function,
we can get separate pileups for each region (as a RaggedArray),
and analyse those pileups using numpy-functionality.
"""

import bionumpy as bnp
import npstructures as nps
import numpy as np

reads = bnp.open("example_data/alignments.bed").read_chunks()
regions = bnp.open("example_data/test2.bed").read()
chrom_sizes = bnp.open("example_data/small.chrom.sizes").read()
chrom_sizes = {str(name): size for name, size in zip(chrom_sizes.name, chrom_sizes.size)}

reads_by_chromosome = bnp.groupby(reads, "chromosome")
regions_by_chromosome = bnp.groupby(regions, "chromosome")

for (chromosome, chromosome_reads), (_, chromosome_regions) in zip(reads_by_chromosome, regions_by_chromosome):
    # Get read pileup, i.e. number of reads per base in this chromosome
    # (we convert the pileup, which is a  to a regular numpy array with to_array())
    read_pileup = bnp.arithmetics.get_pileup(chromosome_reads, chrom_sizes[chromosome]).to_array()

    # Subset the pileup on the regions
    subset = nps.ragged_slice(read_pileup, chromosome_regions.start, chromosome_regions.stop)

    # Subset is now a RaggedArray. Each row contains the pileup for each region
    # We can easily e.g. get the max pileup for each region using axis=-1
    max_per_region = np.max(subset, axis=-1)
    print("")
    print("Chromosome", chromosome)
    print("Regions analysed:", chromosome_regions)
    print("Max pileup within each region:")
    print(max_per_region)
    print("Mean pileup within each region:")
    print(np.mean(subset, axis=-1))