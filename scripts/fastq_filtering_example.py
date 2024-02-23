import numpy as np
import bionumpy as bnp
from bionumpy import streamable

def test(file="example_data/big.fq.gz", out_file="example_data/big_filtered.fq.gz"):
    with bnp.open(out_file, 'w') as out_file:
        for reads in bnp.open(file).read_chunks():
            min_quality_mask = reads.quality.min(axis=-1) > 1
            max_quality_mask = reads.quality.mean(axis=-1) > 10
            mask = min_quality_mask & max_quality_mask
            print(f'Filtering reads: {len(reads)} -> {mask.sum()}')
            out_file.write(reads[mask])


if __name__ == "__main__":
    test()