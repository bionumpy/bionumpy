import dataclasses
from .intervals import sort_intervals
import numpy as np


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

def call_peaks(filename, output_folder):
    reads = bnp.open(filename)
    fragment_pileups = get_fragment_pileup(reads)
    tmp_writer.passthrough_write(fragment_pileups)
    summaries = 
