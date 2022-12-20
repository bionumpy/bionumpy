k = 5
from collections import defaultdict
counts = defaultdict(int)

with open(snakemake.input[0]) as f:
    for i, line in enumerate(f):
        if line.startswith(">"):
            continue

        kmers = [
            line.strip()[i:i+k] for i in range(0, len(line)-k)
        ]
        for kmer in kmers:
            counts[kmer] += 1

        if i % 10000 == 0:
            print(i, "lines processed")


sorted_kmers = sorted(counts)
with open(snakemake.output[0], "w") as f:
    for kmer in sorted_kmers:
        f.write("%s %d\n" % (kmer, counts[kmer]))

