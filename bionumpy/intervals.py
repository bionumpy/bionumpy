import numpy as np
import dataclasses
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
    all_intervals = np.concatenate([intervals_a, intervals_b])
    all_intervals = all_intervals[np.argsort(all_intervals.start, kind="mergesort")]
    ends = np.sort(all_intervals.end, kind="mergesort")
    mask = ends[:-1] > all_intervals.start[1:]
    result = all_intervals[1:][mask]
    result.end = ends[:-1][mask]
    return result

@ChromosomeMap()
def extend(intervals, both=None, forward=None, reverse=None, left=None, right=None):
    directed = (forward is not None) or (reverse is not None)
    undirected = (left is not None) or (right is not None)
    assert sum([both is not None, directed, undirected]) == 1
    if both is not None:
        return dataclasses.replace(intervals, 
                                   start=intervals.start-both,
                                   end=intervals.end+both)
    if undirected:
        starts = interval.start-left if left is not None else interval.start
        ends = interval.end+right if right is not None else interval.end
        dataclasses.replace(intervals, 
                            start=starts,
                            end=ends)
    if directed:
        if forward is None:
            forward = 0
        if reverse is None:
            reverse = 0
        return dataclasses.replace(
            intervals,
            start = np.where(interval.strand=="+",
                             intervals.start-reverse,
                             intervals.start-forward),
            end = np.where(interval.strand=="+",
                             intervals.end+forward,
                             intervals.start+reverse)
            )

def pileup(intervals):
    positions = np.concatenate((intervals.start,
                                 intervals.end))
    args = np.argsort(positions, kind="mergesort")
    values = np.where(args >= len(intervals), -1, 1)
    np.cumsum(values, out=values)
    positions = positions[args]
    mask = positions[1:] != positions[:-1]
    print(positions, values, mask)
    kept_positions = positions[:-1][mask]
    values = values[:-1][mask]
    return BedGraph(intervals.chromosome, 
                    positions, 
                    np.append(kept_positions[1:], positions[-1]),
                    np
                    

# .sort(kind="mergesort", key="start")
#     starts = np.concatenate([intervals_a.start, intervals_b.start])
#     ends = np.concatenate([intervals_a.end, intervals_b.end])
#     starts.sort(kind="mergesort")
#     ends.sort(kind="mergesort")
#     mask = ends[:-1] > starts[1:]
#     return intervals_a.__class__(intervals_a.chromosome
#     return np.maximum(, 0)
# 
