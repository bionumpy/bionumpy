import bionumpy as bnp
from bionumpy.bam import alignment_to_interval

alignments = bnp.open("example_data/test.bam").read()
outfile = bnp.open("example_data/converted_alignments.bed", "w")
outfile.write(alignment_to_interval(alignments))
outfile.close()
