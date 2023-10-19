import bionumpy as bnp
from bionumpy.io.one_line_buffer import TwoLineFastaBuffer


def reverse_complement(input_filename, output_filename):
    with bnp.open(output_filename, "w", buffer_type=TwoLineFastaBuffer) as outfile:
        for chunk in bnp.open(input_filename, buffer_type=TwoLineFastaBuffer).read_chunks():
            outfile.write(bnp.replace(chunk, sequence=bnp.sequence.get_reverse_complement(chunk.sequence)))


def _test():
    reverse_complement('../benchmarks/results/dna_sequences/length150_nreads500000.fa', 'tmp.fa')


if __name__ == '__main__':
    import sys

    reverse_complement(sys.argv[1], sys.argv[2])
