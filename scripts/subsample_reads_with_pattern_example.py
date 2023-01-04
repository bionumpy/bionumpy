import bionumpy as bnp
import numpy as np


def subsample(reads, pattern):
    matches = bnp.sequence.string_matcher.match_string(reads.sequence, pattern)
    return reads[np.sum(matches, axis=1) > 0]


def test():
    in_file = bnp.open("example_data/big.fq.gz")
    out_file = bnp.open("out.fa.gz", "w")
    for chunk in in_file.read_chunks():
        out_file.write(subsample(chunk, "ACT"))
