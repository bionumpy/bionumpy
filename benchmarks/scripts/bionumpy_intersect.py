import bionumpy as bnp
import dataclasses

chrom_sizes = bnp.open(snakemake.input.chrom_sizes).read()
a = bnp.open(snakemake.input.a, buffer_type=bnp.Bed6Buffer).read()
b = bnp.open(snakemake.input.b, buffer_type=bnp.Bed6Buffer).read()


print("Sorting")

chromosome_encoding = bnp.encodings.string_encodings.StringEncoding(chrom_sizes.name)
a = dataclasses.replace(a, chromosome=chromosome_encoding.encode(a.chromosome))
b = dataclasses.replace(b, chromosome=chromosome_encoding.encode(b.chromosome))

# a = bnp.arithmetics.sort_intervals(a, sort_order=chrom_sizes.name.tolist())
# b = bnp.arithmetics.sort_intervals(b, sort_order=chrom_sizes.name.tolist())
print("DOne sorting")

#ms = bnp.streams.MultiStream(chrom_sizes, a=a, b=b)
result = bnp.arithmetics.global_intersect(a, b)
with bnp.open(snakemake.output[0], "w", buffer_type=bnp.Bed6Buffer) as f:
    f.write(result)
