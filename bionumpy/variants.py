def is_snp(variant):
    return (variant.ref_seq.shape.lengths == 1) & (variant.alt_seq.shape.lengths == 1)
