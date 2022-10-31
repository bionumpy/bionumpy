.. _intervals:

===========
 Intervals
===========

Working with intervals is a central part of many bioinformatics analyses. In BioNumPy, an interval is any object (`bnpdataclass`) having a `start` and `stop` attribute and optionally a `strand` attribute. We will typically get intervals by reading them from files such as `.bed`, `.gtf` or `.narrowPeaks`. Or by converting alignments to reference intervals by `alignment.to_interval`. The functionality in BioNumPy for handling intervals can be divided into two operating on intervals as independent objects, or operations with some sort of dependencies:


Independent Operations
======================

Independent operations covers any thing we want to do with a set of intervals that can be done one interval at a time. This covers simple geometric operations such as shifting and resizing; filtering based on properties of intervals. Since both the start and stop (and strand) attribute are numpy arrays, we can implement alot of functionality simply by using numpy indexing and ufuncs:

    >>> import bionumpy as bnp
    >>> from dataclasses import replace
    >>> intervals = bnp.open("example_data/small_interval.bed").read()
    >>> extended_right = replace(intervals, stop=intervals.stop+10)
    >>> shifted = replace(intervals, start=intervals.start+5, stop=intervals.stop+5)
    >>> small = intervals[(intervals.stop-intervals.start)<50]

Even though it is possible to write these kinds of functions using simple numpy functionality, it is often better to wrap them in functions with instructive names. In the `bnp.interval` module, a set of utility funcitons ffor working with intervals is provided.

Operations with dependencies
============================
Some queries cannot be answered just by looking at intervals as independent entries. For instance: "All regions" covered by at least two intervals". In order to answer such questions, we need to look at the intervals as a coherent set. In BioNumPy, such queries are answered using `npstructures.RunLengthArray". 


Genomic Intervals
=================
When working with genomic intervals, we often want to only deal intervals from one chromsome at the time. In order to do this, we can use the `bnp.group_by` function on interval entries.

    >>> import bionumpy as bnp
    >>> intervals = bnp.open("example_data/small_interval.bed").read()
    >>> print(intervals)
    Interval with 50 entries
                   chromosome                    start                     stop
                            0                       13                       18
                            0                       37                       46
                            0                       62                       83
                            0                      105                      126
                            0                      129                      130
                            1                        3                       21
                            1                       41                       65
                            1                       91                      114
                            1                      131                      153
                            1                      157                      168
    >>> for chromosome, data in bnp.groupby(intervals, "chromosome"):
    ...     print(f"---Chromosome: {chromosome}---")
    ...     print(data)
    ---Chromosome: 0---
    Interval with 5 entries
                   chromosome                    start                     stop
                            0                       13                       18
                            0                       37                       46
                            0                       62                       83
                            0                      105                      126
                            0                      129                      130
    ---Chromosome: 1---
    Interval with 10 entries
                   chromosome                    start                     stop
                            1                        3                       21
                            1                       41                       65
                            1                       91                      114
                            1                      131                      153
                            1                      157                      168
                            1                      174                      201
                            1                      213                      230
                            1                      240                      268
                            1                      290                      315
                            1                      319                      339
    ---Chromosome: 2---
    Interval with 15 entries
                   chromosome                    start                     stop
                            2                        2                       16
                            2                       44                       49
                            2                       77                      101
                            2                      108                      127
                            2                      135                      154
                            2                      163                      165
                            2                      173                      177
                            2                      201                      214
                            2                      242                      268
                            2                      292                      320
    ---Chromosome: 3---
    Interval with 20 entries
                   chromosome                    start                     stop
                            3                        7                       34
                            3                       58                       82
                            3                       95                      101
                            3                      130                      138
                            3                      150                      170
                            3                      188                      211
                            3                      234                      261
                            3                      283                      302
                            3                      325                      352
                            3                      353                      362
