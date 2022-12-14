import bionumpy as bnp


rule intersect_bionumpy:
    input:
        a="results/intervals/{b}.bed",
        b="results/intervals/{a}.bed",
        chrom_sizes="../example_data/hg38.chrom.sizes"
    output:
        "results/bionumpy/intersect/{a}-vs-{b}.bed"
    benchmark:
        "benchmarks/intersect/bionumpy/{a}-vs-{b}.txt"
    script:
        "../scripts/bionumpy_intersect.py"


rule intersect_bedtools:
    input:
        left="results/intervals/{b}.bed",
        right="results/intervals/{a}.bed"
    output:
        "results/bedtools/intersect/{a}-vs-{b}.bed"
    benchmark:
        "benchmarks/intersect/bedtools/{a}-vs-{b}.txt"
    log:
        "logs/bedtools/intersect/{a}-vs-{b}.log"
    params:
        extra=""
    wrapper:
        "v1.19.2/bio/bedtools/intersect"


rule intersect_pyranges:
    input:
        a = "results/intervals/{b}.bed",
        b = "results/intervals/{a}.bed",
    output:
        "results/pyranges/intersect/{a}-vs-{b}.bed"
    benchmark:
        "benchmarks/intersect/pyranges/{a}-vs-{b}.txt"
    script:
        "../scripts/pyranges_intersect.py"
