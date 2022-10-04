from functools import reduce
import numpy as np
from npstructures import npdataclass
from npstructures.runlengtharray import RunLengthArray


def sum_largest(stream):
    return reduce(lambda a, b: np.pad(a, max(a.size, b.size))+ np.pad(b, max(a.size, b.size)),
                  stream)


def value_hist(graph):
    weights = graph.end-graph.start
    return np.bincount(graph.value, weights=weights)


@npdataclass
class BedGraph:
    chromosome: str
    start: int
    end: int
    value: int


def get_pileup(intervals, size):
    positions = np.concatenate(([0], intervals.start,
                                intervals.end, [size]))
    args = np.argsort(positions, kind="mergesort")
    values = np.where(args >= (len(intervals)+1), -1, 1)
    values[0] = 0
    np.cumsum(values, out=values)
    positions = positions[args]
    mask = np.flatnonzero(positions[1:] == positions[:-1])
    positions = np.delete(positions, mask)
    values = np.delete(values, mask)
    assert positions[-1] == size, (positions[-10:], size)
    return RunLengthArray(positions, values[:-1])


def memory_efficient_pileup(intervals, size):
    intervals.end.sort(kind="mergesort")
    intervals.start.sort(kind="mergesort")
    end_idxs_in_start = np.searchsorted(intervals.start, intervals.end[:-1], side="right")-1
    start_idxs_in_end = np.searchsorted(intervals.end, intervals.start[1:], side="right")-1
    indices = np.arange(len(intervals))
    end_values = end_idxs_in_start-indices
    start_values = start_idxs_in_end-indices
    corresponding_start = np.searchsorted(intervals.start, end_positions)


    positions = np.concatenate(([0], intervals.start,
                                intervals.end, [size]))
    args = np.argsort(positions, kind="mergesort")
    values = np.where(args >= (len(intervals)+1), -1, 1)
    values[0] = 0
    np.cumsum(values, out=values)
    positions = positions[args]
    mask = np.flatnonzero(positions[1:] == positions[:-1])
    positions = np.delete(positions, mask)
    values = np.delete(values, mask)

    #intervals = np.lib.stride_tricks.sliding_window_view(positions, 2)
    #mask = np.flatnonzero(intervals[:, 0] == intervals[:, 1])
    # intervals = np.delete(intervals, mask, axis=0)
    # values = np.delete(values, mask)
    #mask = np.flatnonzero(values[1:] == values[:-1])
    #values = np.delete(values, mask)
    #starts = np.delete(intervals[:, 0], mask+1)
    #ends = np.delete(intervals[:, 1], mask)
    assert positions[-1] == size, (positions[-10:], size)
    return RunLengthArray(positions, values[:-1])
