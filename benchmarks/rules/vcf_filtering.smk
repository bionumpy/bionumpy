
rule download_phased_vcf:
    output:
        "results/vcfs/{size}_phased.vcf.gz"
    shell:
        "wget -O {output} https://github.com/bionumpy/bionumpy-example-data/raw/master/{wildcards.size}_phased.vcf.gz"


rule bcftools_vcf_filtering:
    input:
        "results/vcfs/{size}_phased.vcf.gz"
    output:
        "results/bcftools/vcf_filtering/{size}.filtered.vcf"
    benchmark:
        repeat("benchmarks/vcf_filtering/bcftools/{size}.txt", config["n_benchmark_repeats"])
    log:
        "log/{size}.view.vcf.log"
    #params:
    #    extra="--min-ac 10"
    #wrapper:
    #    "v1.21.0/bio/bcftools/view"
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --min-ac 10 -o {output} {input}
        """


rule bionumpy_vcf_filtering:
    input:
        "results/vcfs/{size}_phased.vcf.gz"
    output:
        "results/bionumpy/vcf_filtering/{size}.filtered.vcf"
    benchmark:
        repeat("benchmarks/vcf_filtering/bionumpy/{size}.txt", config["n_benchmark_repeats"])
    shell:
        """
        python3 -O ../scripts/vcf_allele_frequency_filtering_example.py {input} {output} 10
        """


# not used, similar to pyvcf
rule vcfpy_vcf_filtering:
    input:
        "results/vcfs/{size}_phased.vcf.gz"
    output:
        "results/vcfpy/vcf_filtering/{size}.filtered.vcf"
    benchmark:
        repeat("benchmarks/vcf_filtering/vcfpy/{size}.txt", config["n_benchmark_repeats"])
    run:
        import vcfpy

        reader = vcfpy.Reader.from_path(input[0])
        writer = vcfpy.Writer.from_path(output[0],reader.header)

        for record in reader:
            genotype_calls = record.calls
            allele_count = 0
            for call in record.calls:
                genotype = call.data.get("CT")
                if genotype == "0|1" or genotype == "1|0":
                    allele_count += 1
                elif genotype == "1|1":
                    allele_count += 2

            if allele_count >= 10:
                writer.write_record(record)


rule pyvcf_vcf_filtering:
    input:
        "results/vcfs/{size}_phased.vcf.gz"
    output:
        "results/pyvcf/vcf_filtering/{size}.filtered.vcf"
    benchmark:
        repeat("benchmarks/vcf_filtering/pyvcf/{size}.txt", config["n_benchmark_repeats"])
    conda:
        "../envs/pyvcf.yml"
    script:
        "../scripts/pyvcf_vcf_filtering.py"
