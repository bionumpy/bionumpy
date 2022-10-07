import numpy as np
import dataclasses
from .bedgraph import BedGraph
from .chromosome_map import ChromosomeMap
from .bnpdataclass import bnpdataclass


@bnpdataclass
class RawInterval:
    start: int
    end: int


@ChromosomeMap()
def sort_intervals(intervals):
    args = np.lexsort((intervals.end, intervals.start))
    return intervals[args]


@ChromosomeMap()
def merge_intervals(intervals, distance=0):
    ends = np.maximum.accumulate(intervals.end)
    if distance > 0:
        ends += distance
    valid_start_mask = intervals.start[1:] > intervals[:-1].end
    start_mask = np.concatenate(([True], valid_start_mask))
    end_mask = np.concatenate((valid_start_mask, [True]))
    new_interval = intervals[start_mask]
    new_interval.end = ends[end_mask]
    if distance > 0:
        new_interval.end -= distance
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
    chroms = np.concatenate([intervals.chromosome, intervals.chromosome])
    positions = np.concatenate((intervals.start,
                                intervals.end))
    args = np.argsort(positions, kind="mergesort")
    values = np.where(args >= len(intervals), -1, 1)
    np.cumsum(values, out=values)
    positions = positions[args]
    intervals = np.lib.stride_tricks.sliding_window_view(positions, 2)
    mask = np.flatnonzero(intervals[:, 0] == intervals[:, 1])
    intervals = np.delete(intervals, mask, axis=0)
    values = np.delete(values, mask)
    mask = np.flatnonzero(values[1:] == values[:-1])
    values = np.delete(values, mask)
    starts = np.delete(intervals[:, 0], mask+1)
    ends = np.delete(intervals[:, 1], mask)
    return BedGraph(chroms[:values.size-1],
                    starts, ends, values[:-1])
