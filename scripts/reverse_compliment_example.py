import bionumpy as bnp
import bionumpy.io.one_line_buffer

input_filename = None
output_filname = None


def reverse_complement(input_filename, output_filename):
    sequence_entries = bnp.open(input_filename, buffer_type=bionumpy.io.one_line_buffer.TwoLineFastaBuffer).read_chunks()
    reversed_entries = bnp.sequence.get_reverse_complement(sequence_entries)
    with bnp.open(output_filename, "w", buffer_type=bionumpy.io.one_line_buffer.TwoLineFastaBuffer) as outfile:
        outfile.write(reversed_entries)

def _test():
    assert False
    reverse_complement('benchmarks/results/dna_sequences/length150_nreads500000.fa', 'tmp.fa')

