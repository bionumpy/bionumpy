"""
Meant to work identically as "seqtk gc" (https://github.com/lh3/seqtk/blob/master/seqtk.c#L1501)
For now, not identical.
"""

import bionumpy as bnp
import numpy as np


def _get_continous_regions(a, min_width=20):
    d = np.diff(a, axis=-1, prepend=0, append=0)
    starts = np.nonzero(d == 1)
    ends = np.nonzero(d == -1)
    assert np.all(ends[0] == starts[0])
    assert len(ends[1]) == len(starts[1])

    lengths = ends[1]-starts[1]
    valid = np.nonzero(lengths >= min_width)

    rows = starts[0][valid]
    column_starts = starts[1][valid]
    column_ends = ends[1][valid]

    return rows, column_starts, column_ends


def get_gc_regions(chunk, min_width=20):
    gc_mask = chunk.sequence == "G" or chunk.sequence == "C"
    reads, starts, ends = _get_continous_regions(gc_mask, min_width=min_width)
    print(reads, starts, ends)


def _test():
    f = bnp.open("example_data/big.fa.gz", buffer_type=bnp.TwoLineFastaBuffer)
    get_gc_regions(f.read())


if __name__ == "__main__":
    _test()