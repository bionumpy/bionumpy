rule jellyfish_count:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/jellyfish/kmer_counts/{filename}.jf"
    log:
        "{filename}.jf.log",
    params:
        kmer_length=5,
        size="1G",
        #extra="--canonical",
    threads: 1
    benchmark:
        repeat("benchmarks/kmer_counts/jellyfish/count-{filename}.txt", config["n_benchmark_repeats"])
    wrapper:
        "v1.19.2-20-g6055e791/bio/jellyfish/count"


rule jellyfish_dump:
    input:
        "results/jellyfish/kmer_counts/{prefix}.jf",
    output:
        "results/jellyfish/kmer_counts/{prefix}.csv",
    log:
        "results/jellyfish/kmer_counts/{prefix}.log",
    params:
        extra="-c -t",
    benchmark:
        repeat("benchmarks/kmer_counts/jellyfish/dump-{prefix}.txt", config["n_benchmark_repeats"])
    wrapper:
        "v1.19.2-20-g6055e791/bio/jellyfish/dump"


rule bionumpy_count:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/bionumpy/kmer_counts/{filename}.csv"
    benchmark:
        repeat("benchmarks/kmer_counts/bionumpy/{filename}.txt", config["n_benchmark_repeats"])
    shell:
        'python ../scripts/kmer_counting_example.py {input} {output}'
        #         "../scripts/bionumpy_count_kmers.py"


rule python_count:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/python/kmer_counts/{filename}.csv"
    benchmark:
        repeat("benchmarks/kmer_counts/python/{filename}.txt", config["n_benchmark_repeats"])
    script:
        "../scripts/python_kmer_counting.py"

rule biopython_count:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/biopython/kmer_counts/{filename}.csv"
    benchmark:
        repeat("benchmarks/kmer_counts/biopython/{filename}.txt", config["n_benchmark_repeats"])
    conda:
        "../envs/biopython.yml"
    shell:
        "python3 scripts/biopython_kmer_counting.py {input} {output}"







