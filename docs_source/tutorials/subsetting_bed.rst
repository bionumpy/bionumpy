Getting read pileup for multiple regions
-------------------------------------------

This is a brief example of how you can you can analyse alignments inside a given set of regions. This example assumes that:

* You have a set of alignments as a bed-file (can be converted from bam using bedtools BamToBed).
* You have a set of regions in another bed-file
* For simplicity, we here assumes both files are sorted by chromosome. If not, `bnp.arithmetics.sort_intervals` may be used.

We assume the alignments-file is quite big, so that we don't want to load all of it into memory at once.
Thus, we read the file by chromsome using `bnp.groupby` and create a pileup for each chromosome (a NumPy-array).
This pileup can then easily be analysed, e.g. by slicing it with the regions we are interested in using `npstructures.raggedslice`.
The result is a RaggedArray, which behaves very much like a NumPy matrix. Full example below:


.. testcode::

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

In our small example data, we only have three regions for chromosome 1 and 2 regions on chromosome 2. The output should be:

.. testoutput::

    Chromosome chr1
    Regions analysed: Interval with 3 entries
                   chromosome                    start                     stop
                         chr1                        5                       10
                         chr1                       15                       20
                         chr1                      100                      110
    Max pileup within each region:
    [3 1 0]
    Mean pileup within each region:
    [2.4 1.  0. ]

    Chromosome chr2
    Regions analysed: Interval with 2 entries
                   chromosome                    start                     stop
                         chr2                        1                       10
                         chr2                       25                       35
    Max pileup within each region:
    [1 3]
    Mean pileup within each region:
    [0.55555556 3.        ]


