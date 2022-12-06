from more_itertools import grouper
import numpy as np
n_to_subsample = 1000000
row = 0
n_entries = sum(1 for line in open(snakemake.input[0]))//2
rows = set(np.random.choice(np.arange(n_entries), n_to_subsample, replace=False))
with open(snakemake.output[0], "w") as f:
    for i, entry in enumerate(grouper(open(snakemake.input[0]), 2)):
        if i in rows:
            f.write(''.join(entry))
