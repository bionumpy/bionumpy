

rule sample_fq:
    output:
        "results/dna_sequences/small.fq.gz"
    shell:
        "cp ../example_data/big.fq.gz {output}"


rule download_encode_fastq_file:
    output:
        fq="results/dna_sequences/{id,ENCFF[A-Z0-9]+}.fq.gz",
    shell:
        """
        wget -O {output.fq} https://www.encodeproject.org/files/{wildcards.id}/@@download/{wildcards.id}.fastq.gz
        """


rule download_encode_bed_file:
    output:
        bed="results/bed_files/{id}.bed.gz",
        sorted_bed="results/bed_files/{id}.sorted.bed"
    shell:
        """
        wget -O {output.bed} https://www.encodeproject.org/files/{wildcards.id}/@@download/{wildcards.id}.bed.gz
        zcat {output.bed} | sort -k1,1 -k2,2n > {output.sorted_bed}
        """



rule download_mapped_reads_bed_file:
    output:
        "results/intervals/ENCFF{id}_mapped_reads_{size}.bed"
    shell:
        "wget -O {output} https://github.com/bionumpy/bionumpy-example-data/raw/master/ENCFF{wildcards.id}_mapped_reads_{wildcards.size}.bed"