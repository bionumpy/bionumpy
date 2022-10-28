from npstructures import ragged_slice
from numpy import array, append
import bionumpy as bnp
from bionumpy.intervals import get_boolean_mask


def analyze_within_chromosome(seq_fn, genes_fn):
    chr1 = bnp.open(seq_fn).read()
    genes = bnp.open(genes_fn).read()
    print("Gene-regions: ", genes)
    print("Fasta: ", chr1)
    gc_inside_and_outside(chr1, genes)

def gc_inside_and_outside(chr1, genes):
    gc_inside = get_gc_content(chr1.sequence, genes.start, genes.stop)

    outside_starts = append(array([0]), genes.stop)
    outside_ends = append(genes.start, array([chr1.sequence.shape.size]))
    gc_outside = get_gc_content(chr1.sequence, outside_starts, outside_ends)


    print(f"GC inside: {gc_inside:.2f}, GC outside: {gc_outside:.2f}")

def get_gc_content(sequence, starts, ends):
    seq2 = ragged_slice(sequence, starts, ends)
    aa_counts = {nn: (seq2 == nn).sum() for nn in "ACTG"}
    total_count = sum(aa_counts.values())
    aa_proportions = {nn: aa_counts[nn] / total_count for nn in aa_counts}
    gc_content = sum([aa_proportions[nn] for nn in "GC"])
    return gc_content

analyze_within_chromosome("example_data/gc_test_onechr.fa", "example_data/gc_bedtest_onechr.bed")
