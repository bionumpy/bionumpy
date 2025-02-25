import bionumpy as bnp


rule unique_intersect_bedtools:
    input:
        left="results/intervals/{a}.bed",
        right="results/bed_files/{b}.bed.gz",
    output:
        "results/bedtools/unique_intersect/{a}-vs-{b}.bed"
    benchmark:
        repeat("benchmarks/unique_intersect/bedtools/{a}-vs-{b}.txt", config["n_benchmark_repeats"])
    params:
        extra="-u"
    wrapper:
        "v1.19.2/bio/bedtools/intersect"

rule unique_intersect_bionumpy:
    input:
        a="results/intervals/{a}.bed",
        b="results/bed_files/{b}.bed.gz",
        chrom_sizes="../example_data/hg38.chrom.sizes"
    output:
        "results/bionumpy/unique_intersect/{a}-vs-{b}.bed"
    benchmark:
        repeat("benchmarks/unique_intersect/bionumpy/{a}-vs-{b}.txt", config["n_benchmark_repeats"])
    shell:
        "python ../scripts/unique_intersect_example.py {input} {output}"
        #        "../scripts/bionumpy_unique_intersect.py"


rule intersect_bionumpy:
    input:
        a="results/intervals/{a}.bed",
        b="results/bed_files/{b}.bed.gz",
        chrom_sizes="../example_data/hg38.chrom.sizes"
    output:
        "results/bionumpy/intersect/{a}-vs-{b}.bed"
    benchmark:
        repeat("benchmarks/intersect/bionumpy/{a}-vs-{b}.txt", config["n_benchmark_repeats"])
    script:
        "../scripts/bionumpy_intersect.py"


rule intersect_bedtools:
    input:
        left="results/intervals/{a}.bed",
        right="results/bed_files/{b}.bed.gz"
    output:
        "results/bedtools/intersect/{a}-vs-{b}.bed"
    benchmark:
        repeat("benchmarks/intersect/bedtools/{a}-vs-{b}.txt", config["n_benchmark_repeats"])
    log:
        "logs/bedtools/intersect/{a}-vs-{b}.log"
    params:
        extra=""
    wrapper:
        "v1.19.2/bio/bedtools/intersect"


rule unique_intersect_pyranges:
    input:
        a = "results/intervals/{a}.bed",
        b = "results/bed_files/{b}.bed.gz",
    output:
        "results/pyranges/unique_intersect/{a}-vs-{b}.bed"
    benchmark:
        repeat("benchmarks/unique_intersect/pyranges/{a}-vs-{b}.txt", config["n_benchmark_repeats"])
    conda:
        "../envs/pyranges.yml"
    shell:
        "time python scripts/pyranges_intersect.py {input} {output}"
    #script:
    #    "../scripts/pyranges_intersect.py"
