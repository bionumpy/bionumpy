

rule download_encode_fastq_file:
    output:
        fq="results/dna_sequences/{id}.fq.gz",
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

