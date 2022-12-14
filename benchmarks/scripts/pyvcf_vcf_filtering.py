import gzip
import vcf

reader = vcf.Reader(filename=snakemake.input[0])
writer = vcf.Writer(open(snakemake.output[0], "w"), reader)

for record in reader:
    allele_count = 0
    for call in record.samples:
        genotype = call["GT"]
        if genotype == "0|1" or genotype == "1|0":
            allele_count += 1
        elif genotype == "1|1":
            allele_count += 2

    if allele_count >= 10:
        writer.write_record(record)
