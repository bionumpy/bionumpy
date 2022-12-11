"""
Note: BioNumPy and seqtk are run in "single-pass" mode, i.e. reading whole file into memory and subsampling.
"""
import time

import bionumpy as bnp
import numpy as np


def n_to_subsample(wildcards):
    n_reads = int(wildcards.filename.split(".")[0].split("nreads")[1])
    return n_reads // 2


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


rule seqtk_subsample:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/seqtk/subsample/{filename}.fa"
    params:
        n=n_to_subsample,
        seed=123
    conda:
        "../envs/seqtk.yml"
    benchmark:
        "benchmarks/subsample/seqtk/{filename}.txt"
    #wrapper:
    # wrapper uses pigz to compress
    #    "v1.20.0/bio/seqtk/subsample/se"
    shell:
        "seqtk sample -s 123 -2 {input} {params.n} > {output}"


