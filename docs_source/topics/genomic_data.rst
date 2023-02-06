Genomic Data
============

Analysis of genomic data is best done with the knowledge of how the genome looks like. At the very least, we should know how the names of each chromosome in the genome, and how long the sequence of each chromosome is. This serves several purposes:

* We (and the computer) becomes aware of which reference genome we are using, which can prevent errors from mismatched reference genomes
* Chromosomes without any data gets handled correctly
* We get a speed increase in the algorithms

For these reasons, we always recomend using the genomic data framework for handling genomic data in BioNumPy, not using raw bed, vcf or bdg files directly from the `bnp.open` interface. The main access point of the genomic data interface is the genome object, which can be created directly from a `dict` of chromosome names and sizes; or, more commonly, from a `.chrom.sizes` or `.fa.fai` file, which are tab separated lines of chromosome names and sizes: 

.. testcode::

    import bionumpy as bnp
    genome = bnp.Genome.from_file("example_data/hg38.chrom.sizes")
    print(genome)

.. testoutput::

                   Chromosome                     Size
                         chr1                248956422
                         chr2                242193529
                         chr3                198295559
                         chr4                190214555
                         chr5                181538259
                         chr6                170805979
                         chr7                159345973
                         chr8                145138636
                         chr9                138394717
                        chr10                133797422

From this genome object we can create genomic data from reading files, or from bnpdataclass objects directly. For instance, we can now read a set of intervals from a `.bed` file and get a `GenomicIntervals` object. Because  `GenomicIntervals` objects knows which reference genome they belong to, they are easier and faster to work with than raw 'bnp.Interval' objects. 

.. testcode::

    intervals = genome.read_intervals('example_data/peaks.narrowPeak')
    print(intervals)
    mask = intervals.get_mask()
    print(mask)

.. testoutput::


    Genomic Intervals on ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', '...']:
    Interval with 3 entries
                   chromosome                    start                     stop
                         chr1                  9356548                  9356648
                         chr1                  9358722                  9358822
                         chr1                  9361082                  9361182
    chr1: [ False False False ... False False False]
    chr2: [ False False False ... False False False]
    chr3: [ False False False ... False False False]
    chr4: [ False False False ... False False False]
    chr5: [ False False False ... False False False]
    chr6: [ False False False ... False False False]
    chr7: [ False False False ... False False False]
    chr8: [ False False False ... False False False]
    chr9: [ False False False ... False False False]
    chr10: [ False False False ... False False False]
    .
    .
    .
