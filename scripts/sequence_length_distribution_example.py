import bionumpy as bnp
import numpy as np


def get_sorted_sequence_length_distribution(file_name: str):
    '''Return the sequence lengths found and their frequency'''
    counts = bnp.bincount(chunk.sequence.lengths for chunk in bnp.open(file_name).read_chunks())
    sizes = np.flatnonzero(counts)
    return sizes, counts[sizes]


def test():
    get_sorted_sequence_length_distribution("example_data/big.fq.gz")


if __name__ == "__main__":
    import sys
    file_name = sys.argv[1]
    out_file = sys.argv[2]
    sizes, counts = get_sorted_sequence_length_distribution(file_name)
    with open(out_file, "w") as f:
        f.write("\n".join(["%d %d" % (size, count) for size, count in zip(sizes, counts)]) + "\n")
