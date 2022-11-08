import bionumpy as bnp
import pyfaidx
from bionumpy.encodings.alphabet_encoding import DNAEncoding
from bionumpy.encoded_array import EncodedArray, EncodedRaggedArray
from bionumpy.datatypes import SequenceEntry, Variant, Interval, VCFEntry
import numpy as np

np.random.seed(4321)


def simulate_genome(chrom_sizes):
    letters = np.random.randint(0, 4, size=sum(chrom_sizes))
    letters = EncodedArray(letters, DNAEncoding)
    sequences = EncodedRaggedArray(letters, chrom_sizes)
    return SequenceEntry([str(i) for i, _ in enumerate(chrom_sizes)], sequences)


def simulate_snps(bases, chrom_name, rate=0.05):
    # bases = sequence_entries.sequence.ravel()
    positions = np.flatnonzero(np.random.choice([0, 1], size=bases.size, p=[1-rate, rate]))
    ref_seq = bases[positions]
    alt_seq = EncodedArray(np.random.randint(0, 4, size=ref_seq.size), DNAEncoding)
    print(ref_seq, alt_seq)
    mask = ref_seq != alt_seq
    print(mask)
    return VCFEntry([chrom_name]*int(mask.sum()), positions[mask], ["."]*int(mask.sum()),
                    EncodedRaggedArray(ref_seq[mask], np.ones(ref_seq[mask].size, dtype=int)),
                    EncodedRaggedArray(alt_seq[mask], np.ones(alt_seq[mask].size, dtype=int)))


def simulate_intervals(chrom_name, chrom_size):
    diffs = np.random.randint(1, 30, size=chrom_size//30)
    positions = np.cumsum(diffs)
    positions = positions[positions < chrom_size]
    if len(positions) % 2 == 1:
        positions = positions[:-1]
    positions = positions.reshape(-1, 2)
    return Interval([chrom_name]*len(positions), positions[:, 0], positions[:, 1])


entries = simulate_genome([300*i for i in range(1, 5)])
snps = np.concatenate(
    [simulate_snps(entry.sequence, entry.name.to_string())
     for entry in entries])
intervals = np.concatenate(
    [simulate_intervals(entry.name.to_string(), len(entry.sequence)) for entry in entries])

bnp.open("example_data/small_genome.fa", "w").write(entries)
pyfaidx.Faidx("example_data/small_genome.fa")
print(snps)
bnp.open("example_data/few_variants.vcf", "w").write(snps)
bnp.open("example_data/small_interval.bed", "w").write(intervals)
