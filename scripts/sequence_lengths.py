import bionumpy as bnp


@bnp.streamable()
def get_sequence_lengths(sequence_entries):
    return sequence_entries.sequence.shape.lengths


def get_line_length_distribution(fastq_file_name):
    sequence_entries = bnp.open(fastq_file_name).read_chunks()
    return bnp.bincount(get_sequence_lengths(sequence_entries))


def test():
    print(get_line_length_distribution("example_data/reads.fq"))

if __name__ == "__main__":
    test()
