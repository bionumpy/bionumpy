import numpy as np
from .chromosome_map import ChromosomeMap


@ChromosomeMap()
def filter_on_intervals(entry, sorted_intervals):
    if len(sorted_intervals)==0:
        mask = np.full(entry.position.shape, False)
    else:
        starts, ends = (sorted_intervals.start, sorted_intervals.end)
        idx = np.searchsorted(starts, entry.position, side="right")-1
        idx = np.minimum(idx, starts.size-1)
        mask = (entry.position >= starts[idx]) & (entry.position < ends[idx])
    return entry[mask]


@ChromosomeMap()
def get_snps(variants):
    snps = variants[variants.is_snp()]
    snps.ref_seq = snps.ref_seq.ravel()
    snps.alt_seq = snps.alt_seq.ravel()
    return snps
