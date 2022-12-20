
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
        "benchmarks/kmer_counts/jellyfish/count-{filename}.txt"
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
        "benchmarks/kmer_counts/jellyfish/dump-{prefix}.txt"
    wrapper:
        "v1.19.2-20-g6055e791/bio/jellyfish/dump"


rule bionumpy_count:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/bionumpy/kmer_counts/{filename}.csv"
    benchmark:
        "benchmarks/kmer_counts/bionumpy/{filename}.txt"
    script:
        "../scripts/bionumpy_count_kmers.py"


rule python_count:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/python/kmer_counts/{filename}.csv"
    benchmark:
        "benchmarks/kmer_counts/python/{filename}.txt"
    script:
        "../scripts/python_kmer_counting.py"







