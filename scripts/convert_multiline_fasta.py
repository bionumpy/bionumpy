import bionumpy as bnp
import bionumpy.io.one_line_buffer


def convert_to_two_line_fasta(input_file_name, output_file_name):
    sequence_entries = bnp.open(input_file_name).read_chunks()
    with bnp.open(output_file_name, "w", buffer_type=bionumpy.io.one_line_buffer.TwoLineFastaBuffer) as out_file:
        out_file.write(sequence_entries)


def test():
    convert_to_two_line_fasta("example_data/small_genome.fa", "example_data/two_line_genome.fa")


if __name__ == "__main__":
    test()
