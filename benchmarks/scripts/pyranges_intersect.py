import pyranges

a = pyranges.read_bed(snakemake.input.a)
b = pyranges.read_bed(snakemake.input.b)

intersection = a.intersect(b)
intersection.to_bed(snakemake.output[0])
