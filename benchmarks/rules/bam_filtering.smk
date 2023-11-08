rule pysam_filter_bam:
    input:
        'results/bams/{size}.bam'
    output:
        'results/pysam/filtered_bams/{size}.filtered.bam'
    run:
        import pysam
        samfile = pysam.AlignmentFile(input[0], "rb")
        pairedreads = pysam.AlignmentFile(output[0], "wb", template=samfile)
        for read in samfile:
            if read.mapq==60:
                pairedreads.write(read)
        
        pairedreads.close()
        samfile.close()

rule bionumpy_filter_bam:
    input:
        'results/bams/{size}.bam'
    output:
        'results/bionumpy/filtered_bams/{size}.filtered.bam'
    run:
        import bionumpy as bnp
        with bnp.open(output[0], 'w') as f:
            for chunk in bnp.open(input[0]).read_chunks():
                f.write(chunk[chunk.mapq==60])
