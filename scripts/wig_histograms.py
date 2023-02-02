import bionumpy as bnp
import numpy as np
from bionumpy.arithmetics import GenomicTrack
f = bnp.open('/home/knut/Downloads/out.wig')
chrom_sizes = bnp.open('/home/knut/Data/hg38.chrom.sizes').read()
chrom_sizes = dict(zip(chrom_sizes.name.tolist(), chrom_sizes.size))
track = GenomicTrack.from_bedgraph(np.concatenate(list(f.read_chunks())), chrom_sizes)
print(track)
