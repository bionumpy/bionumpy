import numpy as np
from ..encodings import DNAEncoding
from ..datatypes import Variant
from ..encoded_array import EncodedArray
from ..encoded_array import as_encoded_array
from ..sequence import get_reverse_complement, count_encoded
from ..streams import streamable
from ..sequence.lookup import Lookup
import logging

logger = logging.getLogger(__name__)


def get_kmer_indexes(position, flank=1):
    return np.add.outer(position, np.arange(-flank, flank + 1))


class SNPEncoding:
    lookup = Lookup(np.full((4, 4), 255, dtype=np.uint8), DNAEncoding)
    lookup["C", "AGT"] = np.arange(3)
    lookup["G", "TCA"] = np.arange(3)
    lookup["T", "ACG"] = np.arange(3)+3
    lookup["A", "TGC"] = np.arange(3)+3
    text = np.array([f"C>{c}" for c in "AGT"] + [f"T>{c}" for c in "ACG"])

    @classmethod
    def to_string(cls, encoded):
        return cls.text[encoded]

    @classmethod
    def encode(cls, snp):
        values = cls.lookup[snp.ref_seq, snp.alt_seq]
        assert not np.any(values==255)
        return EncodedArray(cls.lookup[snp.ref_seq, snp.alt_seq], cls)

    @classmethod
    def decode(cls, encoded):
        pass


class MutationTypeEncoding:
    def __init__(self, flank, encoding=DNAEncoding):
        k = flank*2+1
        self.k = k
        self.h = 4 ** np.arange(k)
        self.h[k // 2 + 1 :] = self.h[k // 2 : -1]
        self.h[k // 2] = 0
        self.h = self.h[::-1]
        self._encoding = encoding

    def encode(self, kmer: EncodedArray, snp: Variant) -> np.ndarray:
        kmer = as_encoded_array(kmer, self._encoding)
        snp.ref_seq = as_encoded_array(snp.ref_seq.ravel(), self._encoding)
        snp.alt_seq = as_encoded_array(snp.alt_seq.ravel(), self._encoding)
        assert not np.any(snp.ref_seq == snp.alt_seq)
        assert kmer.shape[-1] == self.k, (kmer.shape, self.k)
        assert np.all(kmer[..., self.k//2] == snp.ref_seq), (kmer, snp.ref_seq)
        forward_mask = (snp.ref_seq == "C") | (snp.ref_seq == "T")
        kmer = np.where(forward_mask[:, None], kmer, get_reverse_complement(kmer))
        kmer = kmer.raw()
        kmer_hashes = np.dot(kmer, self.h)
        snp_hashes = SNPEncoding.encode(snp).raw()
        return EncodedArray((kmer_hashes + 4 ** (self.k - 1) * snp_hashes),
                            self)

    def decode(self, encoded):
        snp = SNPEncoding.decode(encoded >> (2 * (self.k - 1)))
        chars = (encoded >> (2 * np.arange(self.k - 1))) & 3
        kmer = "".join(chr(b) for b in self._encoding._decode(chars))
        return kmer[: self.k // 2] + "[" + snp + "]" + kmer[self.k // 2 :]

    def to_string(self, encoded):
        snp = SNPEncoding.to_string(encoded >> (2 * (self.k - 1)))
        chars = (encoded >> (2 * np.arange(self.k - 1))) & 3
        kmer = "".join(chr(b) for b in self._encoding._decode(chars))[::-1]
        return kmer[: self.k // 2] + "[" + snp + "]" + kmer[self.k // 2 :]

    def get_labels(self):
        return [self.to_string(c) for c in np.arange(4**(self.k-1)*6)]


def encode_snps(kmer, alt_seq, true_ref_seq=None):
    kmer = as_encoded_array(kmer, DNAEncoding)
    alt_seq = as_encoded_array(alt_seq.ravel(), DNAEncoding)
    k = kmer.shape[-1]
    ref_seq = kmer[..., k//2]
    if true_ref_seq is not None:
        assert np.all(ref_seq.ravel() == true_ref_seq.ravel())
    forward_mask = (ref_seq == 'C') | (ref_seq == 'T')
    kmer = np.where(forward_mask[:, None], kmer, get_reverse_complement(kmer))
    # alt_seq = np.where(forward_mask, alt_seq, get_reverse_complement(kmer))
    snp_code = SNPEncoding.lookup[ref_seq, alt_seq]
    encoding = MutationTypeEncoding(k//2)
    kmer_hashes = np.dot(kmer.raw(), encoding.h)
    return EncodedArray((kmer_hashes + 4 ** (k - 1) * snp_code),
                        encoding)


from ..genomic_data import GenomicSequence, GenomicLocation


def count_mutation_types_genomic(variants: GenomicLocation, reference: GenomicSequence, flank=1):
    snp_mask = (variants.get_data_field('alt_seq').shape[-1] ==1) & (variants.get_data_field('ref_seq').shape[-1]==1)
    snps = variants[snp_mask]
    ref_seq = snps.get_data_field('ref_seq')
    windows = variants.get_windows(flank)
    kmers = reference[windows].to_numpy_array()
    hashes = encode_snps(kmers, snps.get_data_field('alt_seq'), ref_seq)
    return count_encoded(hashes)


@streamable(reduction=sum)
def count_mutation_types(variants, reference, flank=1):
    snps = variants[variants.is_snp()]
    if len(snps) == 0:
        return 0
    snps = snps[np.argsort(snps.position)]
    mnv_mask = (snps.position[1:] == (snps.position[:-1]+1))
    # mask = np.append(mask, False) | np.insert(mask, 0, False)
    # snps = snps[~mask]
    reference = as_encoded_array(reference)
    kmer_indexes = get_kmer_indexes(snps.position, flank=flank)
    kmers = reference[kmer_indexes]
    mask = np.any((kmers=="n") | (kmers=="N"), axis=-1)
    kmers = kmers[~mask]
    snps = snps[~mask]
    hashes = MutationTypeEncoding(flank).encode(kmers, snps)
    if not hasattr(snps, "genotypes"):
        return count_encoded(hashes)
    has_snps = (snps.genotypes>0)
    masks = (has_snps[1:] & has_snps[:-1]) & mnv_mask[:, np.newaxis]
    masks = np.pad(masks, [(1, 0), (0, 0)]) | np.pad(masks, [(0, 1), (0, 0)])
    has_snps &= ~masks
    counts = []
    for genotype in (snps.genotypes.T>0):
        idxs = np.flatnonzero(genotype)
        if len(idxs):
            # assert np.all(snps.position[idxs[:-1]] < snps.position[idxs[1:]]), snps.position[idxs]
            mask = (snps.position[idxs[:-1]]+1) >= (snps.position[idxs[1:]])
            mask = np.append(mask, False) | np.insert(mask, 0, False)
            idxs = idxs[~mask]
        counts.append(count_encoded(hashes[idxs]))
    counts = EncodedCounts.vstack(counts)
    # assert np.all(counts.counts.sum(axis=-1) == (snps.genotypes>0).sum(axis=0))
    # ~((np.append(mask, False) | np.insert(mask, 0, False)
    return counts

