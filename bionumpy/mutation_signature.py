import numpy as np
from .encodings import SimpleEncoding, ACTGEncoding, BaseEncoding
from .chromosome_map import chromosome_map, ChromosomeMap
import logging
logger = logging.getLogger(__name__)

def get_kmer_indexes(position, flank=2):
    relative = np.concatenate((np.arange(-flank, 0),
                               np.arange(1, flank+1)))
    return position + relative


class SNPEncoding:
    lookup = np.zeros((256, 256), dtype=np.uint8)
    lookup[ord("C")][[ord(c) for c in "AGT"]] = np.arange(3)
    lookup[ord("G")][[ord(c) for c in "TCA"]] = np.arange(3)
    lookup[ord("T")][[ord(c) for c in "ACG"]] = 3+np.arange(3)
    lookup[ord("A")][[ord(c) for c in "TGC"]] = 3+np.arange(3)

    text = np.array([f"C>{c}" for c in "AGT"] + [f"T->{c}" for c in "ACG"])

    @classmethod
    def to_string(cls, encoded):
        return cls.text[encoded]
    
    @classmethod
    def from_snp(cls, snp):
        return cls.lookup[snp.ref_seq, snp.alt_seq]

    @classmethod
    def to_snp(cls, encoded):
        pass


class MutationSignatureEncoding:
    @classmethod
    def from_kmers_and_snp(cls, kmer, snp):
        k = kmer.shape[-1]
        h = 4**np.arange(k)
        kmer_hashes = np.sum(ACTGEncoding.from_bytes(kmer)*h, axis=-1)
        snp_hashes = SNPEncoding.from_snp(snp)
        return kmer_hashes + 4**k*snp_hashes

    @classmethod
    def to_string(cls, encoded, k):
        snp = SNPEncoding.to_string(encoded>>(2*k))
        chars = (encoded>>(2*np.arange(k))) & 3

        kmer = "".join(chr(b) for b in ACTGEncoding.to_bytes(chars))
        return kmer[:k//2]+"["+snp+"]"+kmer[k//2:]
        return snp+ ":" +kmer
        kmer_bytes = ACTGEncoding.to_bytes(encoded)

@ChromosomeMap(reduction=sum)
def get_kmers(variants, genotypes, intervals, reference, flank):
    
    assert np.all(reference[variants.position] == variants.ref_seq[:, 0:1].ravel())
    snps = variants[variants.is_snp()]
    snps.ref_seq = snps.ref_seq.ravel()# to_numpy_array().ravel()
    snps.alt_seq = snps.alt_seq.ravel()# to_numpy_array().ravel()
    assert np.all(reference[snps.position] == snps.ref_seq), (reference[snps.position], snps.ref_seq)
    if intervals is not None:
        valid_indexes = np.flatnonzero(intervals.in_intervals(snps.position))
        snps = snps[valid_indexes]
        if genotypes is not None:
            genotypes = genotypes[valid_indexes]
    kmer_indexes = get_kmer_indexes(snps.position[:, None], flank=flank)
    kmers = reference[kmer_indexes]

    forward_mask = (snps.ref_seq == ord("C")) | (snps.ref_seq==ord("T"))
    forward_idxs = np.flatnonzero(forward_mask)
    reverse_idxs = np.flatnonzero(~forward_mask)
    forward_hashes = MutationSignatureEncoding.from_kmers_and_snp(kmers[forward_idxs], snps[forward_idxs])
    reverse_hashes = MutationSignatureEncoding.from_kmers_and_snp(
        BaseEncoding.complement(kmers[reverse_idxs, ::-1]), snps[reverse_idxs])
    all_hashes = np.zeros_like(snps.position)
    all_hashes[forward_mask] = forward_hashes
    all_hashes[~forward_mask] = reverse_hashes

    n_hashes = 4**(flank*2)*6

    if genotypes is not None:
        count_matrix = np.array([
            np.bincount(all_hashes, weights=genotypes[:, sample] > 0, minlength=n_hashes)
            for sample in range(genotypes.shape[-1])], dtype=int)
    else:
        count_matrix = np.bincount(all_hashes, minlength=n_hashes)
    return count_matrix
            
"""
def test_get_kmers():
    seq = np.array([ord(c) for c in "AAAGCAAAATGCAAATTCAAAGAA"],dtype=np.uint8)
    intervals = SortedIntervals(np.array([[0, 100]], dtype=int))
    snps = SNP(np.ones((4, 1)), 4+6*np.arange(4), np.array([ord(c) for c in "CGTA"], dtype=np.uint8),
               np.array([ord(c) for c in "TAGC"], dtype=np.uint8))
    print(get_kmers((snps, np.ones((4, 1))), intervals, seq, 1))
"""
