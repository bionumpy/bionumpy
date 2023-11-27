from collections import Counter
from Bio import SeqIO


def count_kmers(input_filename, output_filename, k=5):
    counts = Counter(kmer for dna_record in SeqIO.parse(input_filename, 'fasta') for kmer in get_kmers(k, dna_record.seq))
    with open(output_filename, "w") as f:
        f.writelines(f'{kmer} {counts[kmer]}\n' for kmer in sorted(counts))


def get_kmers(k, seq):
    return (
        seq[pos:pos + k] for pos in range(0, len(seq) - k+1)
    )


def test():
    i = '../results/dna_sequences/length150_nreads200.fa'
    o = 'tmp.csv'
    count_kmers(i, o)


if __name__ == '__main__':
    import sys
    #test()
    i, o = sys.argv[1:3]
    count_kmers(i, o)
