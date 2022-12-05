from .intervals import get_boolean_mask
from ..streams import MultiStream, streamable
import numpy as np


@streamable(sum)
def get_contingency_table(intervals_a, intervals_b, sequence_length):
    boolean_a = get_boolean_mask(intervals_a, sequence_length)
    boolean_b = get_boolean_mask(intervals_b, sequence_length)
    return np.array([[np.sum(boolean_a & boolean_b), np.sum(boolean_a & ~boolean_b)],
                     [np.sum(~boolean_a & boolean_b), np.sum(~boolean_a & ~boolean_b)]])


def forbes(chromosome_sizes, intervals_a, intervals_b):
    ms = MultiStream(chromosome_sizes, a=intervals_a, b=intervals_b)
    ((a, b), (c, d)) = get_contingency_table(ms.a, ms.b, ms.lengths)
    N = (a+b+c+d)
    return a*N/((a+b)*(a+c))


def jaccard(chromosome_sizes, intervals_a, intervals_b):
    ms = MultiStream(chromosome_sizes, a=intervals_a, b=intervals_b)
    ((a, b), (c, d)) = get_contingency_table(ms.a, ms.b, ms.lengths)
    N = (a+b+c+d)
    return a/(N-d)
