import numpy as np
import cupy as cp

import bionumpy as bnp
from bionumpy.io.file_buffers import TwoLineFastaBuffer

bnp.set_backend(cp)

fasta_filename = "../../data/fa/testreads20m.fa"
chunk_size = 5000000

chunk_generator = bnp.open(fasta_filename, buffer_type=TwoLineFastaBuffer).read_chunks(min_chunk_size=chunk_size)
for chunk, _ in zip(chunk_generator, range(5)):
    encoded_sequence = bnp.as_encoded_array(chunk.sequence, bnp.encodings.alphabet_encoding.ACTGEncoding)
    kmers = bnp.kmers.fast_hash(encoded_sequence, 31, bnp.encodings.alphabet_encoding.ACTGEncoding).ravel()
    print(type(kmers))
    print(kmers)

