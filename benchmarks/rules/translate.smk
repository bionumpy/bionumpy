import bionumpy as bnp


rule translate_bionumpy:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/bionumpy/translate/{filename}.fa"
    benchmark:
        "benchmarks/translate/bionumpy/{filename}.txt"
    script:
        "../scripts/bionumpy_translate.py"


rule translate_biopython:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biopython/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/biopython/{name}.txt"
    conda:
        "../envs/biopython.yml"
    script:
        "../scripts/biopython_translate.py"


rule translate_biostrings:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biostrings/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/biostrings/{name}.txt"
    script:
        "scripts/reverse_complement_biostrings.R"


rule translate_python:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/python/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/python/{name}.txt"
    script:
        "../scripts/python_translate.py"


rule translate_biotite:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biotite/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/biotite/{name}.txt"
    conda:
        "../envs/biotite.yml"
    script:
        "../scripts/biotite_translate.py"
