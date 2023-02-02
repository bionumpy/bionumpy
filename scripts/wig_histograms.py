import bionumpy as bnp
import numpy as np
import matplotlib.pyplot as plt
from bionumpy.arithmetics import GenomicTrack
f = bnp.open('/home/knut/Downloads/out.wig')
chrom_sizes = bnp.open('/home/knut/Data/hg38.chrom.sizes').read()
chrom_sizes = dict(zip(chrom_sizes.name.tolist(), chrom_sizes.size))
track = GenomicTrack.from_bedgraph(f.read_chunks(), chrom_sizes)
values, edges = np.histogram(track, bins=100, range=(0, 3)).compute()
print(values, edges)

plt.bar(edges[:-1], values); plt.show()
