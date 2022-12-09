from bionumpy.sequence import get_kmers, count_encoded
from bionumpy.streams import streamable
from bionumpy.io.dump_csv import dump_csv
import bionumpy as bnp


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

result = Result(kmers.alphabet, kmers.counts)
print(result)
output_stream.write(result)
