from functools import reduce
from npstructures import npdataclass


def sum_largest(stream):
    return reduce(lambda a, b: np.pad(a, max(a.size, b.size))+ np.pad(b, max(a.size, b.size)),
                  stream)
    
# @bnp_broadcast(sum_largest)
def value_hist(graph):
    weights = graph.end-graph.start
    return np.bincount(graph.value, weights=weights)

@npdataclass
class BedGraph:
    chromosome: str
    start: int
    end: int
    value: int
