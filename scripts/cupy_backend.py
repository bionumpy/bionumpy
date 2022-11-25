import cupy as cp
import sys
import bionumpy as bnp
import bionumpy.encoded_array_functions
from bionumpy.io.file_buffers import TwoLineFastaBuffer

bnp.set_backend(cp)


def run(fasta_filename="example_data/small.fa"):
    chunk_size = 500000
    chunk_generator = bnp.open(fasta_filename, buffer_type=TwoLineFastaBuffer).read_chunks(min_chunk_size=chunk_size)
    for chunk, _ in zip(chunk_generator, range(5)):
        encoded_sequence = bionumpy.encoded_array_functions.as_encoded_array(chunk.sequence, bnp.encodings.alphabet_encoding.ACTGEncoding)
        kmers = bionumpy.sequence.kmers.get_kmers(encoded_sequence, 31).ravel()
        print(type(kmers))
        print(kmers)


if __name__ == "__main__":
    run(sys.argv[1])
