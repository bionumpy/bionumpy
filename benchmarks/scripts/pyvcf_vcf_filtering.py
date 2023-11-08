import gzip
import vcf

reader = vcf.Reader(filename=snakemake.input[0])
writer = vcf.Writer(open(snakemake.output[0], "w"), reader)

for record in reader:
    if record.INFO['AC'][0] >= 10:
        writer.write_record(record)
