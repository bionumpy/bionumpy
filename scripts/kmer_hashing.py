import sys
import bionumpy as bnp

for chunk in bnp.open(sys.argv[1]):
    kmers = bnp.sequence.get_kmers(chunk.sequence, 31)
    print(kmers.ravel())
