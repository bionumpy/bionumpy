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


@ChromosomeMap(reduction=sum)
def count_overlap(intervals_a, intervals_b):
    starts = np.concatenate([intervals_a.start, intervals_b.start])
    ends = np.concatenate([intervals_a.end, intervals_b.end])
    starts.sort(kind="mergesort")
    ends.sort(kind="mergesort")
    return np.sum(np.maximum(ends[:-1]-starts[1:], 0))


@ChromosomeMap()
def intersect(intervals_a, intervals_b):
    starts = np.concatenate([intervals_a.start, intervals_b.start])
    ends = np.concatenate([intervals_a.end, intervals_b.end])
    starts.sort(kind="mergesort")
    ends.sort(kind="mergesort")
    return np.maximum(ends[:-1]-starts[1:], 0)

