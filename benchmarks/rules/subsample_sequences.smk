"""
Note: BioNumPy and seqtk are run in "single-pass" mode, i.e. reading whole file into memory and subsampling.
"""
import time

import bionumpy as bnp
import numpy as np

n_to_subsample = 100

rule bionumpy_subsample:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/bionumpy/subsample/{filename}.fa"
    benchmark:
        "benchmarks/subsample/bionumpy/{filename}.txt"
    run:
        t = time.perf_counter()
        reads = bnp.open(input[0]).read()
        rows = np.random.choice(np.arange(len(reads.sequence)), n_to_subsample, replace=False)
        subsample = reads[rows]
        bnp.open(output[0], "w").write(subsample)
        print(time.perf_counter()-t)


rule python_subsample:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/python/subsample/{filename}.fa"
    benchmark:
        "benchmarks/subsample/python/{filename}.txt"
    run:
        row = 0

        lines = open(input[0]).readlines()
        rows = np.random.choice(np.arange(len(lines)//2), n_to_subsample, replace=False)

        subsample = [lines[row]+lines[row+1] for row in rows]
        with open(output[0], "w") as f:
            f.write(''.join(subsample))


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


