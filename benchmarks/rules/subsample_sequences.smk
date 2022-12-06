"""
Note: BioNumPy and seqtk are run in "single-pass" mode, i.e. reading whole file into memory and subsampling.
"""
import time

import bionumpy as bnp
import numpy as np

n_to_subsample = 1000000


rule bionumpy_subsample:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/bionumpy/subsample/{filename}.fa"
    benchmark:
        "benchmarks/subsample/bionumpy/{filename}.txt"
    script:
        "../scripts/bionumpy_subsample.py"


rule python_subsample:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/python/subsample/{filename}.fa"
    benchmark:
        "benchmarks/subsample/python/{filename}.txt"
    script:
        "../scripts/python_subsample.py"
#         from more_itertools import grouper
#         row = 0
#         n_entries = sum(1 for line in open(input[0]))
#         rows = set(np.random.choice(np.arange(n_entries), n_to_subsample, replace=False))
#         with open(output[0], "w") as f:
#             for i, entry in enumerate(grouper(open(input[0]), 2)):
#                 if i in rows:
#                     f.write(''.join(entry))

rule seqtk_subsample:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/seqtk/subsample/{filename}.fa"
    log:
        "results/seqtk/subsample/{filename}.log"
    params:
        n=n_to_subsample,
        seed=123
    benchmark:
        "benchmarks/subsample/seqtk/{filename}.txt"
    wrapper:
        "v1.20.0/bio/seqtk/subsample/se"


