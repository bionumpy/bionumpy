import numpy as np
from .encodings import DNAEncoding, BaseEncoding
from .variants import is_snp
from .datatypes import Variant
from .encoded_array import as_encoded_array, EncodedArray
from .dna import reverse_compliment
from .chromosome_map import ChromosomeMap
from .counter import count_encoded, EncodedCounts
from .lookup import Lookup
import logging

logger = logging.getLogger(__name__)




def get_kmer_indexes(position, flank=2):
    return position[..., np.newaxis] + np.arange(-flank, flank + 1)


class SNPEncoding:
    lookup = Lookup(np.zeros((4, 4), dtype=np.uint8), DNAEncoding)
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
        return cls.lookup[snp.ref_seq, snp.alt_seq]

    @classmethod
    def decode(cls, encoded):
        pass


class MutationTypeEncoding:
    def __init__(self, k, encoding=DNAEncoding):
        self.k = k
        self.h = 4 ** np.arange(k)
        self.h[k // 2 + 1 :] = self.h[k // 2 : -1]
        self.h[k // 2] = 0
        self._encoding = encoding

    def encode(self, kmer: EncodedArray, snp: Variant) -> np.ndarray:
        kmer = as_encoded_array(kmer, self._encoding)
        snp.ref_seq = as_encoded_array(snp.ref_seq.ravel(), self._encoding)
        snp.alt_seq = as_encoded_array(snp.alt_seq.ravel(), self._encoding)
        assert kmer.shape[-1] == self.k, (kmer.shape, self.k)
        assert np.all(kmer[..., self.k//2] == snp.ref_seq), (kmer, snp.ref_seq)
        forward_mask = (snp.ref_seq == "C") | (snp.ref_seq == "T")
        kmer = np.where(forward_mask[:, None], kmer, reverse_compliment(kmer))
        kmer = kmer.raw()
        kmer_hashes = np.dot(kmer, self.h)
        snp_hashes = SNPEncoding.encode(snp)
        
        return (kmer_hashes + 4 ** (self.k - 1) * snp_hashes)

    def decode(self, encoded):
        snp = SNPEncoding.decode(encoded >> (2 * (self.k - 1)))
        chars = (encoded >> (2 * np.arange(self.k - 1))) & 3
        kmer = "".join(chr(b) for b in self._encoding.decode(chars))
        return kmer[: self.k // 2] + "[" + snp + "]" + kmer[self.k // 2 :]

    def to_string(self, encoded):
        snp = SNPEncoding.to_string(encoded >> (2 * (self.k - 1)))
        chars = (encoded >> (2 * np.arange(self.k - 1))) & 3
        kmer = "".join(chr(b) for b in self._encoding.decode(chars))
        return kmer[: self.k // 2] + "[" + snp + "]" + kmer[self.k // 2 :]

    def get_labels(self):
        return [self.to_string(c) for c in np.arange(4**(self.k-1)*6)]


@ChromosomeMap(reduction=sum)
def count_mutation_types(variants, reference, flank=1):
    snps = variants[is_snp(variants)]
    reference = as_encoded_array(reference, DNAEncoding)
    # snps.ref_seq = as_encoded_array(snps.ref_seq.ravel(), DNAEncoding)
    # snps.alt_seq = as_encoded_array(snps.alt_seq.ravel(), DNAEncoding)
    kmer_indexes = get_kmer_indexes(snps.position, flank=flank)
    print(kmer_indexes)
    kmers = reference[kmer_indexes]
    # forward_mask = (snps.ref_seq == "C") | (snps.ref_seq == "T")
    # kmers = EncodedArray(
    #     np.where(
    #         forward_mask[:, None], kmers, reverse_compliment(kmers)
    #     ), DNAEncoding)
    signature_encoding = MutationTypeEncoding(flank * 2 + 1)
    all_hashes = signature_encoding.encode(kmers, snps)
    all_hashes = EncodedArray(all_hashes, encoding=signature_encoding)
    n_hashes = 4 ** (flank * 2) * 6
    if not hasattr(snps, "genotypes"):
        return count_encoded(all_hashes)
    return EncodedCounts.concatenate([count_encoded(all_hashes, weights=genotype>0)
                                      for genotype in snps.genotypes.T])
