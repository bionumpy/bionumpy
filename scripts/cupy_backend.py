import cupy as cp
import sys
import bionumpy as bnp
bnp.set_backend(cp)
import bionumpy.encoded_array
from bionumpy.io.file_buffers import TwoLineFastaBuffer



def run(fasta_filename="../example_data/big.fq.gz"):
    chunk_size = 5000
    buffer = TwoLineFastaBuffer
    buffer = bnp.FastQBuffer
    chunk_generator = bnp.open(fasta_filename, buffer_type=buffer).read_chunks(min_chunk_size=chunk_size)
    for chunk, _ in zip(chunk_generator, range(5)):
        #encoded_sequence = bionumpy.encoded_array.as_encoded_array(chunk.sequence, bnp.encodings.alphabet_encoding.ACTGEncoding)
        s = bnp.change_encoding(chunk.sequence, bnp.DNAEncoding)
        kmers = bionumpy.sequence.kmers.get_kmers(s, 31).ravel()
        print(type(kmers))
        print(kmers)


if __name__ == "__main__":
    run()
