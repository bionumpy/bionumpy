import bionumpy as bnp
from bionumpy.multiline_buffer import MultiLineFastaBuffer
import sys
MultiLineFastaBuffer.n_characters_per_line = 80
f = bnp.open(sys.argv[2], "w")
f.write(bnp.open(sys.argv[1]))
f.close()
