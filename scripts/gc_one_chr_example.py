from typing import Tuple

import bionumpy as bnp
from bionumpy.arithmetics import get_boolean_mask


def analyze_within_chromosome(seq_filename: str, genes_filename: str):
    chr1 = bnp.open(seq_filename).read()[0]
    genes = bnp.open(genes_filename).read()
    print("Gene-regions: ", genes)
    print("Fasta: ", chr1)
    gc_inside, gc_outside = gc_inside_and_outside(chr1.sequence, genes)
    print(f"GC inside: {gc_inside:.2f}, GC outside: {gc_outside:.2f}")
    return gc_inside, gc_outside


def gc_inside_and_outside(chr1_sequence: bnp.EncodedArray, genes: bnp.Interval)-> Tuple[int, int]:
    gc_inside = get_gc_content(chr1_sequence, genes)
    gc_outside = get_gc_content(chr1_sequence, ~get_boolean_mask(genes, len(chr1_sequence)))
    return gc_inside, gc_outside


def get_gc_content(sequence, intervals):
    selected_seq = sequence[intervals]
    nn_counts = {nn: (selected_seq == nn).sum() for nn in "ACTG"}
    gc_count = sum([nn_counts[nn] for nn in "GC"])
    return gc_count / sum(nn_counts.values())

if __name__ == '__main__':
    result = analyze_within_chromosome("example_data/gc_test_onechr.fa", "example_data/gc_bedtest_onechr.bed")
    assert result == (0.6, 0.2)

