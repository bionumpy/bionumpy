import bionumpy as bnp

fastq_entries_stream = bnp.open("example_data/reads.fq")
for fastq_entries in fastq_entries_stream:
    print(fastq_entries)

variants_stream = bnp.open("example_data/variants.vcf")
for chromosome, variants in variants_stream:
    print(variants)
