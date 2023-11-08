rule pysam_filter_bam:
    input:
        'results/bams/{size}.bam'
    output:
        'results/pysam/filtered_bams/{size}.filtered.bam'
    benchmark:
        "benchmarks/bam_filtering/pysam/{size}.txt"
    conda:
        "../envs/pysam.yml"
    script:
        '../scripts/pysam_filtering.py'

rule bionumpy_filter_bam:
    input:
        'results/bams/{size}.bam'
    output:
        'results/bionumpy/filtered_bams/{size}.filtered.bam'
    benchmark:
        "benchmarks/bam_filtering/bionumpy/{size}.txt"
    run:
        import bionumpy as bnp
        with bnp.open(output[0], 'w') as f:
            for chunk in bnp.open(input[0]).read_chunks():
                f.write(chunk[chunk.mapq==60])
