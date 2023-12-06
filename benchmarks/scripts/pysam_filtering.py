import pysam
samfile = pysam.AlignmentFile(snakemake.input[0], "rb")
pairedreads = pysam.AlignmentFile(snakemake.output[0], "wb", template=samfile)
for read in samfile:
    if read.mapq == 60:
        pairedreads.write(read)

pairedreads.close()
samfile.close()
