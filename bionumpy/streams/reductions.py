from functools import reduce
import numpy as np
from . import streamable


def bincount_reduce(bincount_a, bincount_b):
    if bincount_a.size >= bincount_b.size:
        bincount_a[:bincount_b.size] += bincount_b
        return bincount_a
    bincount_b[:bincount_a.size] += bincount_a
    return bincount_b


bincount = streamable(lambda x: reduce(bincount_reduce, x))(np.bincount)


def histogram_reduce(histograms):
    hist, edge = next(histograms)
    hist = sum((h[0] for h in histograms))+hist
    return hist, edge


histogram = streamable(histogram_reduce)(np.histogram)


@streamable(sum)
def sum_and_n(array, axis=None):
    if axis is None:
        n = array.size
    elif axis == 0:
        n = len(array)
    return np.append(np.sum(array, axis=axis), n)


@streamable()
def _rowmean(array, axis=None):
    return np.mean(array, axis=axis)


def mean(array, axis=None):
    """Streamable version of numpy.mean

    Calculates the sum and count of every chunk and sums them up to
    calculate the mean

    Parameters
    ----------
    array : ArrayLike
        Either ArrayLike or stream of array-likes
    axis : int
        The axis to do mean over

    """
    if (axis is not None) and axis != 0:
        return _rowmean(array, axis)
    t = sum_and_n(array, axis=axis)
    return t[:-1]/t[-1]


def quantile(array, quantiles, axis=None):
    """
    """
    hist = bincount(array)
    cumulative = np.cumsum(hist)
    idxs = np.searchsorted(cumulative, quantiles*cumulative[-1])
    return idxs
