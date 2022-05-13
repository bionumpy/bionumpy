import numpy as np
from .chromosome_map import ChromosomeMap


@ChromosomeMap()
def sort_intervals(intervals):
    args = np.lexsort((intervals.end, intervals.start))
    return intervals[args]


@ChromosomeMap()
def merge_intervals(intervals):
    ends = np.maximum.accumulate(intervals.end)
    valid_start_mask = intervals.start[1:] > intervals[:-1].end
    start_mask = np.concatenate(([True], valid_start_mask))
    end_mask = np.concatenate((valid_start_mask, [True]))
    new_interval = intervals[start_mask]
    new_interval.end = ends[end_mask]
    return new_interval
