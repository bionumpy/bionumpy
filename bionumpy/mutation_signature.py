import numpy as np
from .encodings import SimpleEncoding, ACTGEncoding, BaseEncoding
from .chromosome_map import ChromosomeMap
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

def filter_snps(snps, intervals):
    valid_indexes = np.flatnonzero(intervals.in_intervals(snps.position))
    return snps[valid_indexes]
    
def get_snps(variants):
    snps = variants[variants.is_snp()]
    snps.ref_seq = snps.ref_seq.ravel()# to_numpy_array().ravel()
    snps.alt_seq = snps.alt_seq.ravel()# to_numpy_array().ravel()
    return snps

@ChromosomeMap(reduction=sum)
def get_kmers(variants, reference, flank):
    snps = get_snps(variants)
    assert np.all(reference[snps.position] == snps.ref_seq)
    kmer_indexes = get_kmer_indexes(snps.position[:, None], flank=flank)
    kmers = reference[kmer_indexes]
    forward_mask = (snps.ref_seq == ord("C")) | (snps.ref_seq==ord("T"))
    kmers = np.where(forward_mask[:, None],
                     kmers,
                     BaseEncoding.complement(kmers[:, ::-1]))
    all_hashes = MutationSignatureEncoding.from_kmers_and_snp(kmers, snps)
    n_hashes = 4**(flank*2)*6
    if not hasattr(snps, "genotypes"):
        return np.bincount(all_hashes, minlength=n_hashes)
    return np.array([
        np.bincount(all_hashes, weights=genotypes[:, sample] > 0, minlength=n_hashes)
        for sample in range(genotypes.shape[-1])], dtype=int)
