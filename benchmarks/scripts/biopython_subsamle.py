import numpy as np
import sys
from Bio import SeqIO


def subsample(input_filename, output_filename):
    n_entries = count_entries(input_filename)
    n_to_subsample = n_entries // 2
    rows = set(np.random.choice(np.arange(n_entries), n_to_subsample, replace=False))

    with open(output_filename, "w") as f:
        for i, dna_record in enumerate(SeqIO.parse(input_filename, 'fasta')):
            if i in rows:
                SeqIO.write(dna_record, f, 'fasta-2line')


def count_entries(input_filename):
    return sum(1 for _ in open(input_filename)) // 2


def test():
    i = '../results/dna_sequences/length150_nreads200.fa'
    o = 'tmp.fa'
    subsample(i, o)
    assert count_entries(o) == count_entries(i) // 2


if __name__ == '__main__':
    i, o = sys.argv[1:3]
    subsample(i, o)
# try:
#     i, o = snakemake.input, snakemake.output
# except NameError:
#     pass
# else:
#     subsample(i, o)
# with open(snakemake.output[0], 'w') as aa_fa:
#     for dna_record in SeqIO.parse(snakemake.input[0], 'fasta'):
#         new_record = SeqRecord(
#             dna_record.seq.reverse_complement(),
#             id=dna_record.id, description="")
#     SeqIO.write(new_record, aa_fa, 'fasta-2line')
