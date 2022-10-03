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
    intervals = np.lib.stride_tricks.sliding_window_view(positions, 2)
    mask = np.flatnonzero(intervals[:, 0] == intervals[:, 1])
    intervals = np.delete(intervals, mask, axis=0)
    values = np.delete(values, mask)
    mask = np.flatnonzero(values[1:] == values[:-1])
    values = np.delete(values, mask)
    starts = np.delete(intervals[:, 0], mask+1)
    ends = np.delete(intervals[:, 1], mask)
    return RunLengthArray(np.append(starts, ends[-1]), values[:-1])
