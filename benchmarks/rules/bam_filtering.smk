rule pysam_filter_bam:
    input:
        'results/bams/{size}.bam'
    output:
        'results/pysam/filtered_bams/{size}.filtered.bam'
    benchmark:
        repeat("benchmarks/bam_filtering/pysam/{size}.txt", config["n_benchmark_repeats"])
    conda:
        "../envs/pysam.yml"
    script:
        '../scripts/pysam_filtering.py'

rule samtools_filter_bam:
    input:
        'results/bams/{size}.bam'
    output:
        'results/samtools/filtered_bams/{size}.filtered.bam'
    benchmark:
        repeat("benchmarks/bam_filtering/samtools/{size}.txt", config["n_benchmark_repeats"])
    conda:
        "../envs/samtools.yml"
    shell:
        'samtools view -b -q 60 {input} > {output}'

        
rule bionumpy_filter_bam:
    input:
        'results/bams/{size}.bam'
    output:
        'results/bionumpy/filtered_bams/{size}.filtered.bam'
    benchmark:
        repeat("benchmarks/bam_filtering/bionumpy/{size}.txt", config["n_benchmark_repeats"])
    run:
        import bionumpy as bnp
        with bnp.open(output[0], 'w') as f:
            for chunk in bnp.open(input[0]).read_chunks():
                f.write(chunk[chunk.mapq == 60])
