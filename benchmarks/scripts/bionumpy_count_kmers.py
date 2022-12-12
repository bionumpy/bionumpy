from bionumpy.sequence import get_kmers, count_encoded
from bionumpy.streams import streamable
import bionumpy as bnp
import numpy as np


@streamable(sum)
def count_kmers(sequence_entries):
    sequence = bnp.change_encoding(sequence_entries.sequence, bnp.DNAEncoding)
    kmers = get_kmers(sequence, k=5)
    return count_encoded(kmers, axis=None)


from bionumpy.io.delimited_buffers import DelimitedBuffer


@bnp.bnpdataclass.bnpdataclass
class Result:
    kmer: str
    count: int


class ResultBuffer(DelimitedBuffer):
    dataclass = Result


stream = bnp.open(snakemake.input[0]).read_chunks()
output_stream = bnp.open(snakemake.output[0], "wb", buffer_type=ResultBuffer)
kmers = count_kmers(stream)
alphabet = np.array(kmers.alphabet)
# sort kmers, so that output is identical to jellyfish
sorting = np.argsort(alphabet)
print(alphabet, sorting, kmers.counts[sorting])
result = Result(list(alphabet[sorting]), kmers.counts[sorting])
print(result)
output_stream.write(result)
