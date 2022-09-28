import dataclasses
from bionumpy.intervals import sort_intervals
from bionumpy.bedgraph import value_hist
from scipy.stats import poisson
import bionumpy as bnp
import numpy as np

# Extend the intervals to be size 150. Starts should


def extend_to_size(intervals, fragment_length):
    is_forward = intervals.strand=="+"
    return dataclasses.replace(
        intervals,
        start=np.where(is_forward, intervals.start, intervals.end-fragment_length),
        end=np.where(is_forward, intervals.end, intervals.start+fragment_length))


@chromosome_map()
def get_fragment_pileup(reads, fragment_length):
    fragments = extend_to_size(reads, fragment_length)
    fragments = sort_intervals(fragments, kind="mergesort")
    return pileup(fragments)


def call_peaks(filename, tmp_file, fragment_length, background):
    reads = bnp.open(filename)
    tmp_writer = bnp.open(tmp_file, "w")
    fragment_pileup = get_fragment_pileup(reads)
    fragment_pileup = tmp_writer.passthrough_write(fragment_pileup)
    hist = value_hist(fragment_pileup)
