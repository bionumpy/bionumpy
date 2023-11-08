import gzip
from collections import defaultdict, Counter
from Bio import SeqIO


def sequence_length_distribution(input_filename, output_filename):
    with gzip.open(input_filename, 'rt') as f:
        counts = Counter(len(dna_record.seq) for dna_record in SeqIO.parse(f, 'fastq'))

    with open(output_filename, "w") as f:
        f.writelines(f'{c} {counts[c]}\n' for c in sorted(counts))


def test():
    i = '../results/dna_sequences/small.fq.gz'
    o = 'tmp.fa'
    sequence_length_distribution(i, o)


if __name__ == '__main__':
    # test()
    import sys
    sequence_length_distribution(*sys.argv[1:3])
