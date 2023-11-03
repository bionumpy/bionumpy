import numpy as np
from npstructures import RaggedArray
from warnings import warn
from ..genomic_data import GenomicSequence, GenomicLocation
from ..encodings import DNAEncoding
from ..encoded_array import Encoding
from ..datatypes import Variant
from ..encoded_array import EncodedArray
from ..encoded_array import as_encoded_array
from ..sequence import get_reverse_complement, count_encoded
from ..streams import streamable
from ..sequence.lookup import Lookup
import logging

logger = logging.getLogger(__name__)


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

def encode_snps(kmer, alt_seq, true_ref_seq=None):
    kmer = as_encoded_array(kmer, DNAEncoding)
    if isinstance(kmer, RaggedArray):
        kmer = kmer.to_numpy_array()
    alt_seq = as_encoded_array(alt_seq.ravel(), DNAEncoding)
    k = kmer.shape[-1]
    ref_seq = kmer[..., k//2]
    if true_ref_seq is not None:
        assert np.all(ref_seq.ravel() == as_encoded_array(true_ref_seq, DNAEncoding).ravel())
    forward_mask = (ref_seq == 'C') | (ref_seq == 'T')
    kmer = np.where(forward_mask[:, None], kmer, get_reverse_complement(kmer))
    snp_code = SNPEncoding.lookup[ref_seq, alt_seq]
    encoding = MutationTypeEncoding(k//2)
    kmer_hashes = np.dot(kmer.raw(), encoding.h)
    return EncodedArray((kmer_hashes + 4 ** (k - 1) * snp_code),
                        encoding)

class MutationTypeEncoding(Encoding):

    def __init__(self, flank, encoding=DNAEncoding):
        k = flank*2+1
        self.k = k
        self.h = 4 ** np.arange(k)
        self.h[k // 2 + 1 :] = self.h[k // 2 : -1]
        self.h[k // 2] = 0
        self.h = self.h[::-1]
        self._encoding = encoding
        self.flank = flank

    def encode(self, seq) -> np.ndarray:
        l = seq.shape[-1]
        assert seq.shape[-1] == self.k+4, (seq.shape, seq)
        kmer_idxs = np.concatenate((np.arange(self.flank), [self.flank+1],
                                    np.arange(l-self.flank, l)))
        kmer = seq[..., kmer_idxs]
        kmer = as_encoded_array(kmer, self._encoding)
        ref_seq = kmer[..., self.k//2]
        alt_seq = as_encoded_array(seq[..., self.flank+3], self._encoding)
        kmer_hashes = np.dot(kmer.raw(), self.h)
        snp_hashes = SNPEncoding.lookup[ref_seq, alt_seq]
        return EncodedArray((kmer_hashes + 4 ** (self.k - 1) * snp_hashes),
                            self)

    def from_flanked_snp(self, kmer, alt_seq, ref_seq=None):
        return encode_snps(kmer, alt_seq, ref_seq)
                         
    def __decode(self, encoded):
        snp = SNPEncoding.decode(encoded >> (2 * (self.k - 1)))
        chars = (encoded >> (2 * np.arange(self.k - 1))) & 3
        kmer = "".join(chr(b) for b in self._encoding._decode(chars))
        return kmer[: self.k // 2] + "[" + snp + "]" + kmer[self.k // 2 :]

    def to_string(self, encoded):
        snp = SNPEncoding.to_string(encoded >> (2 * (self.k - 1)))
        chars = (encoded >> (2 * np.arange(self.k - 1))) & 3
        kmer = "".join(chr(b) for b in self._encoding._decode(chars))[::-1]
        return kmer[: self.k // 2] + "[" + snp + "]" + kmer[self.k // 2 :]

    decode = to_string

    def get_labels(self):
        return [self.to_string(c) for c in np.arange(4**(self.k-1)*6)]




def count_mutation_types_genomic(variants: GenomicLocation, reference: GenomicSequence, flank=1, genotyped=False, genotypes=None):
    snp_mask = (variants.get_data_field('alt_seq').shape[-1] == 1) & (variants.get_data_field('ref_seq').shape[-1]==1)
    snps = variants[snp_mask]
    ref_seq = snps.get_data_field('ref_seq')
    windows = snps.get_windows(flank)
    kmers = reference[windows].to_numpy_array()
    mask = ~np.any(kmers == 'N', axis=-1)
    hashes = encode_snps(kmers[mask], snps[mask].get_data_field('alt_seq'), ref_seq[mask])
    if not genotyped and genotypes is None:
        return count_encoded(hashes)
    if genotypes is None:
        genotypes = (snps[mask].get_data_field('genotypes').raw()>0).T
    else:
        genotypes = genotypes[snp_mask][mask].T
    return count_encoded(hashes, genotypes, axis=-1)
