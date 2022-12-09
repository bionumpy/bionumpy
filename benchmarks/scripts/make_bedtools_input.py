f = open(snakemake.output[0], "w")
f.write(snakemake.params.bedtools_commands)
f.close()