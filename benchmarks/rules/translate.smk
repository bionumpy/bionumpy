import bionumpy as bnp


rule translate_bionumpy:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/bionumpy/translate/{filename}.fa"
    benchmark:
        repeat("benchmarks/translate/bionumpy/{filename}.txt", config["n_benchmark_repeats"])
    shell:
        "python ../scripts/translate_example.py {input} {output}"


rule translate_biopython:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biopython/translate/{name}.fa"
    benchmark:
        repeat("benchmarks/translate/biopython/{name}.txt", config["n_benchmark_repeats"])
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
        repeat("benchmarks/translate/biostrings/{name}.txt", config["n_benchmark_repeats"])
    script:
        "scripts/reverse_complement_biostrings.R"


rule translate_python:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/python/translate/{name}.fa"
    benchmark:
        repeat("benchmarks/translate/python/{name}.txt", config["n_benchmark_repeats"])
    script:
        "../scripts/python_translate.py"


rule translate_biotite:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biotite/translate/{name}.fa"
    benchmark:
        repeat("benchmarks/translate/biotite/{name}.txt", config["n_benchmark_repeats"])
    conda:
        "../envs/biotite.yml"
    script:
        "../scripts/biotite_translate.py"
