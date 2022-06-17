import sys
import bionumpy as bnp
from bionumpy.kmers import fast_hash

for chunk in bnp.open(sys.argv[1]):
    kmers = fast_hash(chunk.sequence, 31)
    print(kmers.ravel())
