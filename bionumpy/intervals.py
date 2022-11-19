import dataclasses

import numpy as np
from npstructures import RunLength2dArray, RunLengthArray

from .bedgraph import BedGraph
from .streams.grouped import chromosome_map
from .datatypes import Interval
from .bnpdataclass import bnpdataclass


@bnpdataclass
class RawInterval:
    start: int
    stop: int


def get_pileup(intervals: Interval, chromosome_size: int) -> RunLengthArray:
    """Get the number of intervals that overlap each position of the chromosome/contig

    This uses run length encoded arrays to handle the sparse data that
    we get from intervals.

    Parameters
    ----------
    intervals : Interval
        Intervals on the same chromosome/contig
    chromosome_size : int
        size of the chromsome/contig

    Examples
    --------
    >>> intervals = Interval(["chr1", "chr1", "chr1"], [3, 5, 10], [8, 7, 12])
    >>> pileup = get_pileup(intervals, 20)
    >>> print(pileup)
    [0 0 0 1 1 2 2 1 0 0 1 1 0 0 0 0 0 0 0 0]

    """
    rla = RunLength2dArray.from_intervals(intervals.start, intervals.stop, chromosome_size)
    return rla.sum(axis=0)


def get_boolean_mask(intervals: Interval, chromosome_size: int):
    """Get a boolean mask representing where any inteval hits

    Uses run length encoded binary arrays to represent the areas
    covered by any interval. The mask that is returned supports numpy ufuncs, 
    so that you can run logical operations on them s.a. `& | ~` and also 
    numpy indexing so you can use it to filter positions and intervals.

    Parameters
    ----------
    intervals : Interval
        Intervals on the same chromosome/contig
    chromosome_size : int
        The size of the chromosome/contig

    Examples
    --------

    >>> intervals = Interval(["chr1", "chr1", "chr1"], [3, 5, 10], [8, 7, 12])
    >>> print(intervals)
    Interval with 3 entries
                   chromosome                    start                     stop
                         chr1                        3                        8
                         chr1                        5                        7
                         chr1                       10                       12
    >>> mask = get_boolean_mask(intervals, 20)
    >>> print(mask.astype(int))
    [0 0 0 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 0 0]

    Get complement of the mask:

    >>> complement = ~mask
    >>> print(complement.astype(int))
    [1 1 1 0 0 0 0 0 1 1 0 0 1 1 1 1 1 1 1 1]

    Get the intersections (`&`) and union (`|`) of the mask and another mask

    >>> other_mask = get_boolean_mask(Interval(["chr1"], [9], [15]), 20)
    >>> intersection = mask & other_mask
    >>> print(intersection.astype(int))
    [0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0]
    >>> union = mask | other_mask
    >>> print(union.astype(int))
    [0 0 0 1 1 1 1 1 0 1 1 1 1 1 1 0 0 0 0 0]

    Find wether some positions overlap the mask:
    >>> print(other_mask[intervals.start])
    [False False  True]
    """
    rla = RunLength2dArray.from_intervals(intervals.start, intervals.stop, chromosome_size)
    return rla.any(axis=0)


@chromosome_map()
def sort_intervals(intervals):
    args = np.lexsort((intervals.stop, intervals.start))
    return intervals[args]


@chromosome_map()
def merge_intervals(intervals: Interval, distance: int = 0) -> Interval:
    """Merge a set of sorted intervals

    Merges any overlapping intervals into a single interval.
    If `distance` is specified, will merge any intervals that are closer
    than `distance` apart.

    Parameters
    ----------
    intervals : Interval
        The sorted intervals to be merged
    distance : int
        Max distance for merging

    Return
    --------
    Merged Intervals

    """

    assert np.all(intervals.start[:-1] <= intervals.start[1:]), "merge_intervals requires intervals sorted on start position"
    stops = np.maximum.accumulate(intervals.stop)
    if distance > 0:
        stops += distance
    valid_start_mask = intervals.start[1:] > intervals[:-1].stop
    start_mask = np.concatenate(([True], valid_start_mask))
    stop_mask = np.concatenate((valid_start_mask, [True]))
    new_interval = intervals[start_mask]
    new_interval.stop = stops[stop_mask]
    if distance > 0:
        new_interval.stop -= distance
    return new_interval


@chromosome_map(reduction=sum)
def count_overlap(intervals_a, intervals_b):
    starts = np.concatenate([intervals_a.start, intervals_b.start])
    stops = np.concatenate([intervals_a.stop, intervals_b.stop])
    starts.sort(kind="mergesort")
    stops.sort(kind="mergesort")
    return np.sum(np.maximum(stops[:-1]-starts[1:], 0))


@chromosome_map()
def intersect(intervals_a, intervals_b):
    all_intervals = np.concatenate([intervals_a, intervals_b])
    all_intervals = all_intervals[np.argsort(all_intervals.start, kind="mergesort")]
    stops = np.sort(all_intervals.stop, kind="mergesort")
    mask = stops[:-1] > all_intervals.start[1:]
    result = all_intervals[1:][mask]
    result.stop = stops[:-1][mask]
    return result

@chromosome_map()
def extend(intervals, both=None, forward=None, reverse=None, left=None, right=None):
    directed = (forward is not None) or (reverse is not None)
    undirected = (left is not None) or (right is not None)
    assert sum([both is not None, directed, undirected]) == 1
    if both is not None:
        return dataclasses.replace(intervals, 
                                   start=intervals.start-both,
                                   stop=intervals.stop+both)
    if undirected:
        starts = interval.start-left if left is not None else interval.start
        stops = interval.stop+right if right is not None else interval.stop
        dataclasses.replace(intervals, 
                            start=starts,
                            stop=stops)
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
            stop = np.where(interval.strand=="+",
                             intervals.stop+forward,
                             intervals.start+reverse)
            )


def pileup(intervals):
    chroms = np.concatenate([intervals.chromosome, intervals.chromosome])

    positions = np.concatenate((intervals.start,
                                intervals.stop))
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
    stops = np.delete(intervals[:, 1], mask)
    return BedGraph(chroms[:values.size-1],
                    starts, stops, values[:-1])
