from .intervals import get_boolean_mask
from ..streams import MultiStream, streamable, BnpStream
import numpy as np
from ..datatypes import ChromosomeSize, Interval


@streamable(sum)
def get_contingency_table(intervals_a, intervals_b, sequence_length):
    boolean_a = get_boolean_mask(intervals_a, sequence_length)
    not_a = ~boolean_a
    boolean_b = get_boolean_mask(intervals_b, sequence_length)
    not_b = ~boolean_b
    return np.array([[np.sum(boolean_a & boolean_b), np.sum(boolean_a & not_b)],
                     [np.sum(not_a & boolean_b), np.sum(not_a & not_b)]])


def forbes(chromosome_sizes: ChromosomeSize, intervals_a: Interval, intervals_b: Interval) -> float:
    """Computes the Forbes similarity index for two sets of intervals.

    Parameters
    ----------
    chromosome_sizes : ChromosomeSize
        A ChromosomeSizes, typically from reading a chromosome.sizes file with bnp.open()
    intervals_a : Interval
        Must be sorted. Can be read using bnp.open on a bed-file.
    intervals_a : Interval
        Must be sorted. Can be read using bnp.open on a bed-file.

    Returns
    -------
    float
        The forbes similarity index.

    Examples
    --------
    >>> from bionumpy.arithmetics import forbes, sort_intervals
    >>> from bionumpy.datatypes import Interval
    >>> a = Interval.from_entry_tuples([("chr1", 10, 20), ("chr2", 20, 30)])
    >>> b = Interval.from_entry_tuples([("chr2", 15, 25), ("chr1", 10, 40)])
    >>> a_sorted = sort_intervals(a, sort_order=["chr1", "chr2"])
    >>> b_sorted = sort_intervals(b, sort_order=["chr1", "chr2"])
    >>> forbes({"chr1": 100, "chr2": 200}, a_sorted, b_sorted)
    5.625
    """
    ms = MultiStream(chromosome_sizes, a=intervals_a, b=intervals_b)
    ((a, b), (c, d)) = get_contingency_table(ms.a, ms.b, ms.lengths)
    N = (a+b+c+d)
    return a*N/((a+b)*(a+c))


def jaccard(chromosome_sizes: ChromosomeSize, intervals_a: Interval, intervals_b: Interval) -> float:
    """Computes the Jaccard similarity index for two sets of intervals.

    Parameters
    ----------
    chromosome_sizes : ChromosomeSize
        A ChromosomeSizes, typically from reading a chromosome.sizes file with bnp.open()
    intervals_a : Interval
        Must be sorted. Can be read using bnp.open on a bed-file.
    intervals_a : Interval
        Must be sorted. Can be read using bnp.open on a bed-file.

    Returns
    -------
    float
        The forbes similarity index.

    Examples
    --------
    See forbes for examples.
    """
    ms = MultiStream(chromosome_sizes, a=intervals_a, b=intervals_b)
    ((a, b), (c, d)) = get_contingency_table(ms.a, ms.b, ms.lengths)
    N = (a+b+c+d)
    return a/(N-d)
