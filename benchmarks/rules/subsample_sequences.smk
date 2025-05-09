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
        repeat("benchmarks/subsample/bionumpy/{filename}.txt", config["n_benchmark_repeats"])
    shell:
        'python3 ../scripts/subsample_reads.py {input} {output}'
        #"../scripts/bionumpy_subsample.py"


rule python_subsample:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/python/subsample/{filename}.fa"
    benchmark:
        repeat("benchmarks/subsample/python/{filename}.txt", config["n_benchmark_repeats"])
    script:
        "../scripts/python_subsample.py"

rule biopython_subsample:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/biopython/subsample/{filename}.fa"
    benchmark:
        repeat("benchmarks/subsample/biopython/{filename}.txt", config["n_benchmark_repeats"])
    conda:
        "../envs/biopython.yml"
    shell:
        'python3 scripts/biopython_subsamle.py {input} {output}'

#         "../scripts/biopython_subsample.py"


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
        repeat("benchmarks/subsample/seqtk/{filename}.txt", config["n_benchmark_repeats"])
    #wrapper:
    # wrapper uses pigz to compress
    #    "v1.20.0/bio/seqtk/subsample/se"
    shell:
        "seqtk sample -s 123 -2 {input} {params.n} > {output}"


