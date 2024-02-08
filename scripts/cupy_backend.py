import cupy as cp
import sys
import bionumpy as bnp
import bionumpy.io.fastq_buffer
from bionumpy.io.parser import NumpyFileReader

bnp.set_backend(cp)
import bionumpy.encoded_array
from bionumpy.io import TwoLineFastaBuffer
from isal import igzip


def run(fasta_filename="../example_data/big.fq.gz"):
    chunk_size = 5000
    buffer = TwoLineFastaBuffer
    buffer = bionumpy.io.fastq_buffer.FastQBuffer
    chunk_generator = bnp.open(fasta_filename, buffer_type=buffer).read_chunks(min_chunk_size=chunk_size)
    for chunk, _ in zip(chunk_generator, range(5)):
        #encoded_sequence = bionumpy.encoded_array.as_encoded_array(chunk.sequence, bnp.encodings.alphabet_encoding.ACTGEncoding)
        d = chunk.get_data()
        print(d)
        continue
        s = bnp.change_encoding(chunk.sequence, bnp.DNAEncoding)
        kmers = bionumpy.sequence.kmers.get_kmers(s, 31).ravel()
        print(type(kmers))
        print(kmers)


def run2(filename="../example_data/big.fq.gz"):
    open_func = igzip.open
    buffer = bionumpy.io.fastq_buffer.FastQBuffer
    reader = NumpyFileReader(open_func(filename, "rb"), buffer)
    reader.set_prepend_mode()  # important for performance

    for chunk in reader:
        print(type(chunk.get_data().sequence.raw()))


def run3(filename="../example_data/big.fq.gz"):
    buffer = bionumpy.io.fastq_buffer.FastQBuffer
    reader = bnp.open(filename, buffer_type=buffer)

    for chunk in reader:
        print(chunk.sequence.raw())


if __name__ == "__main__":
    run3()
