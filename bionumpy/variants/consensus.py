import numpy as np

from bionumpy import replace, EncodedArray, VCFEntry


def apply_variants_to_sequence(sequence: EncodedArray, variants: VCFEntry) -> EncodedArray:
    """
    Apply variants to a sequence by replacing the reference sequence with the alternative sequence.
    Works only for variants where ref sequence and alt sequence have the same length.

    Parameters
    ----------
    sequence : EncodedArray
        The sequence to apply the variants to
    variants : VCFEntry
        The variants to apply

    Returns
    -------
    EncodedArray
        The sequence with the variants applied

    """
    seq = sequence.copy()
    assert np.all(seq[variants.position] == variants.ref_seq.ravel()), (seq[variants.position], seq[variants.position+1], seq[variants.position-1],
                                                                        variants.ref_seq.ravel(), variants.position)
    seq[variants.position] = variants.alt_seq.ravel()
    return seq


def apply_variants(sequence_entries, variants):
    """
    Wrapper around `apply_variants_to_sequence` that applies variants to multiple sequences.
    """
    assert np.all(variants.alt_seq.lengths == 1)
    return replace(sequence_entries, sequence=[
        apply_variants_to_sequence(
            entry.sequence,
            variants[variants.chromosome==entry.name]
        ) for entry in sequence_entries])

