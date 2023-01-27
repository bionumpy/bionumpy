import bionumpy as bnp
from bionumpy.arithmetics.geometry import Geometry


chrom_sizes = bnp.open("example_data/hg38.chrom.sizes").read()
geometry = Geometry.from_chrom_sizes(chrom_sizes)
print(geometry)
print(repr(geometry))

tf = bnp.open("example_data/ctcf.bed.gz").read()
tf2 = bnp.open("example_data/znf263.bed.gz").read()

# jaccard
print(geometry.jaccard(tf, tf2))

reads = bnp.open("example_data/many_alignments.bam", buffer_type=bnp.io.bam.BamIntervalBuffer).read()
pileup = geometry.get_pileup(reads)
intervals_pileup = pileup.get_intervals(tf)

