from bionumpy import TwoLineFastaBuffer
from bionumpy.sequence import get_kmers, count_encoded
import bionumpy as bnp
from bionumpy.io.delimited_buffers import DelimitedBuffer


def count_kmers(sequence_entries):
    sequence = bnp.as_encoded_array(sequence_entries.sequence, bnp.DNAEncoding)
    kmers = get_kmers(sequence, k=5)
    return count_encoded(kmers, axis=None)


def count_file(input_file):
    stream = bnp.open(input_file, buffer_type=TwoLineFastaBuffer).read_chunks()
    return sum(count_kmers(chunk) for chunk in stream)


def main(input_file, output_file):
    kmers = count_file(input_file)
    with open(output_file, 'w') as f:
        f.writelines(f'{kmer}\t{count}\n'
                     for kmer, count
                     in sorted(zip(kmers.alphabet, kmers.counts)))


def _test_profiling():
    main('../benchmarks/results/dna_sequences/length150_nreads10000000.fa', 'tmp.csv')
    assert open('tmp.csv').readline().strip() == 'AAAAA	1426301'


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])