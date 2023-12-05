import numpy as np

from bionumpy import replace


def apply_variants_to_sequence(sequence, variants):
    seq = sequence.copy()
    seq[variants.position] = variants.alt_seq.ravel()
    return seq

def apply_variants(sequence_entries, variants):
    assert np.all(variants.alt_seq.lengths == 1)
    return replace(sequence_entries, sequence=[apply_variants_to_sequence(entry.sequence, variants[variants.chromosome==entry.name]) for entry in sequence_entries])

