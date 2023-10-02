import pyranges
a = pyranges.read_bed(snakemake.input.a)
b = pyranges.read_bed(snakemake.input.b)
j = a.stats.jaccard(b)
with open(snakemake.output[0], "w") as f:
    f.write(str(j) + "\n")
