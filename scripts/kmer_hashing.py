import sys
import bionumpy as bnp
from bionumpy.kmers import fast_hash
from bionumpy.encodings import ACTGEncoding
for chunk in bnp.open(sys.argv[1]):
    kmers = fast_hash(chunk.sequence, 31, ACTGEncoding)
    print(kmers.ravel())
