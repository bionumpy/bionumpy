import dataclasses
from bionumpy.bam import alignment_to_interval
from bionumpy.groupby import groupby
from bionumpy.intervals import sort_intervals
from bionumpy.bedgraph import value_hist, get_pileup
from bionumpy.delimited_buffers import Bed6Buffer
from scipy.stats import poisson
from scipy.special import pdtrc
from bionumpy.parser import chunk_lines
import matplotlib.pyplot as plt
import bionumpy as bnp
import numpy as np
import logging

logging.basicConfig(level=logging.INFO)
# Extend the intervals to be size 150. Starts should


def extend_to_size(intervals, fragment_length, size):
    is_forward = intervals.strand == "+"
    start = np.where(is_forward, intervals.start, np.maximum(intervals.end-fragment_length, 0))
    end = np.where(is_forward, intervals.end, np.minimum(intervals.start+fragment_length, size))
    return dataclasses.replace(
        intervals,
        start=start,
        end=end)


def get_fragment_pileup(reads, fragment_length, size):
    logging.info("Getting fragment pileup")
    fragments = extend_to_size(reads, fragment_length, size)
    print(fragments)
    fragments = fragments[np.argsort(fragments.start, kind="mergesort")]
    return get_pileup(fragments, size)


def get_control_pileup(reads, window_sizes, n_reads, size):
    mid_points = (reads.start+reads.end)//2
    pileup = n_reads/size
    for window_size in window_sizes:
        start = np.maximum(mid_points-window_size//2, 0)
        end = np.minimum(mid_points+window_size//2, size)
        windows = dataclasses.replace(reads, start=start, end=end)
        new_pileup= get_pileup(windows, size)
        print(np.mean(np.diff(new_pileup.values)))
        print(new_pileup.values[-1])
        pileup = np.maximum(new_pileup, pileup/window_size)
        print(np.mean(np.diff(pileup.values)))
    return pileup


def logsf(count, mu):
    return np.log(pdtrc(count, mu))


def get_p_values(intervals, fragment_length, chrom_size):
    intervals.strand = intervals.strand.ravel()
    fragment_pileup = get_fragment_pileup(intervals, fragment_length, chrom_size)
    control = fragment_length*get_control_pileup(intervals, [1000, 10000], len(intervals), chrom_size)
    p_values = logsf(fragment_pileup, control)
    plt.plot(p_values.starts, p_values.values); plt.show()


intervals = bnp.open("/home/knut/Data/ENCFF296OGN.bed", buffer_type=Bed6Buffer).read_chunks()
# intervals = alignment_to_interval(reads)
# first_chunk = next(chunk_lines(iter(intervals), 100000))
grouped = groupby(intervals, "chromosome")
chrom, first_chunk = next(iter(grouped))
chrom_size = 248956422
get_p_values(first_chunk, 300, chrom_size)
