import numpy as np

import bionumpy as bnp


def test_read_fasta(fasta_file, genome_sequence):
    genome = bnp.open(fasta_file).read()
    bnp.util.testing.assert_bnpdataclass_equal(genome_sequence, genome)


def test_read_bed(bed_file, intervals):
    intervals = bnp.open(bed_file).read()
    bnp.util.testing.assert_bnpdataclass_equal(intervals, intervals)


def test_fastq(fastq_file, reads):
    reads = bnp.open(fastq_file).read()
    bnp.util.testing.assert_bnpdataclass_equal(reads, reads)


def test_interval_intersection(bed_file, bed_file2, genome):
    mask = genome.read_intervals(bed_file).get_mask()
    mask2 = genome.read_intervals(bed_file2).get_mask()
    intersection = mask & mask2
    union = mask | mask2
    xor = (mask & ~mask2) | (~mask & mask2)
    assert xor.sum() == union.sum() - intersection.sum()


def test_pileup(genome, bed_file2):
    pileup = genome.read_intervals(bed_file2).get_pileup()
    intervals = bnp.open(bed_file2).read()
    assert np.sum(pileup) == np.sum(intervals.stop - intervals.start)


def test_read_bedgraph(bedgraph_file, pileup):
    bdg = bnp.open(bedgraph_file).read()
    bdg_sum = np.sum(np.abs(bdg.stop - bdg.start) * bdg.value)
    assert pileup.sum() == bdg_sum
