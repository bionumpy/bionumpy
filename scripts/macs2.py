import dataclasses
from bionumpy.bam import alignment_to_interval
from bionumpy.groupby import groupby
from bionumpy.intervals import sort_intervals
from bionumpy.bedgraph import value_hist, get_pileup
from scipy.stats import poisson
import matplotlib.pyplot as plt
import bionumpy as bnp
import numpy as np

# Extend the intervals to be size 150. Starts should


def extend_to_size(intervals, fragment_length, size):
    is_forward = intervals.strand == "+"
    print(is_forward.shape, intervals.start.shape, intervals.end.shape)
    start = np.where(is_forward, intervals.start, np.maximum(intervals.end-fragment_length, 0))
    print(start)
    end = np.where(is_forward, intervals.end, np.minimum(intervals.start+fragment_length, size))
    print(end)
    return dataclasses.replace(
        intervals,
        start=start,
        end=end)


def get_fragment_pileup(reads, fragment_length, size):
    fragments = extend_to_size(reads, fragment_length, size)
    print(fragments)
    fragments = fragments[np.argsort(fragments.start)]
    return get_pileup(fragments, size)


def get_control_pileup(reads, window_sizes, n_reads, size):
    mid_points = (reads.start+reads.end)//2
    pileup = n_reads/size
    for window_size in window_sizes:
        start = np.maximum(mid_points-window_size//2, 0)
        end = np.maximum(mid_points+window_size//2, size)
        windows = dataclasses.replace(reads, start=start, end=end)
        pileup = np.maximum(pileup, get_pileup(windows, size)/window_size)
    return pileup


def call_peaks(treatment_filename, control_filename, fragment_length):
    reads = bnp.open(treatment_filename).read_chunks()
    fragment_pileup = get_fragment_pileup(reads, fragment_length)
    hist = value_hist(fragment_pileup)


all_reads = bnp.open("/home/knut/Data/ENCFF296OGN.bam").read_chunks()
grouped = groupby(all_reads, "chromosome")
chrom, reads = next(iter(grouped))
# chrom, reads = next(all_reads)
intervals = alignment_to_interval(reads)
intervals.strand = intervals.strand.ravel()
print(intervals)
fragment_pileup = get_fragment_pileup(intervals, 300, 15072421)
plt.plot(fragment_pileup.starts, fragment_pileup.values)
plt.show()
control = get_control_pileup(intervals, [1000, 10000], len(intervals), 15072421)
plt.plot(control.starts, control.values);plt.show()
