import numpy as np
from .encodings import ACTGEncoding
from .sequences import as_encoded_sequence_array, as_sequence_array, EncodedArray
from .dna import complement, reverse_compliment
from .chromosome_map import ChromosomeMap
from .counter import count_encoded, EncodedCounts
from . import DNAArray 
import logging

logger = logging.getLogger(__name__)


def get_kmer_indexes(position, flank=2):
    return position[..., np.newaxis] + np.arange(-flank, flank + 1)


class SNPEncoding:
    lookup = np.zeros((4, 4), dtype=np.uint8)
    lookup[int(as_encoded_sequence_array("C", DNAArray).ravel())][as_encoded_sequence_array("AGT", DNAArray)] = np.arange(3)
    lookup[int(as_encoded_sequence_array("G", DNAArray).ravel())][as_encoded_sequence_array("TCA", DNAArray)] = np.arange(3)
    lookup[int(as_encoded_sequence_array("T", DNAArray).ravel())][as_encoded_sequence_array("ACG", DNAArray)] = 3 + np.arange(3)
    lookup[int(as_encoded_sequence_array("A", DNAArray).ravel())][as_encoded_sequence_array("TGC", DNAArray)] = 3 + np.arange(3)

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
    def __init__(self, k, encoding=DNAArray):
        self.k = k
        self.h = 4 ** np.arange(k)
        self.h[k // 2 + 1 :] = self.h[k // 2 : -1]
        self.h[k // 2] = 0
        self._encoding = encoding

    def encode(self, kmer, snp):
        kmer = as_encoded_sequence_array(kmer, self._encoding)
        assert kmer.shape[-1] == self.k, (kmer.shape, self.k)
        kmer = np.asarray(kmer)
        kmer_hashes = np.dot(kmer, self.h)
        snp_hashes = SNPEncoding.encode(snp)
        return (kmer_hashes + 4 ** (self.k - 1) * snp_hashes)

    def decode(self, encoded):
        snp = SNPEncoding.decode(encoded >> (2 * (self.k - 1)))
        chars = (encoded >> (2 * np.arange(self.k - 1))) & 3
        kmer = "".join(chr(b) for b in self._encoding.encoding.decode(chars))
        return kmer[: self.k // 2] + "[" + snp + "]" + kmer[self.k // 2 :]

    def to_string(self, encoded):
        snp = SNPEncoding.to_string(encoded >> (2 * (self.k - 1)))
        chars = (encoded >> (2 * np.arange(self.k - 1))) & 3
        kmer = "".join(chr(b) for b in self._encoding.encoding.decode(chars))
        return kmer[: self.k // 2] + "[" + snp + "]" + kmer[self.k // 2 :]

    def get_labels(self):
        return [self.to_string(c) for c in np.arange(4**(self.k-1)*6)]
    
def get_encoded_array_class(_encoding):
    class MutationTypeArray(EncodedArray):
        encoding=_encoding
    f = _encoding.encode
    _encoding.encode = lambda *args, **kwargs: f(*args, **kwargs).view(MutationTypeArray)

    return MutationTypeArray

@ChromosomeMap(reduction=sum)
def count_mutation_types(snps, reference, flank=1):
    reference = as_encoded_sequence_array(reference, DNAArray)
    snps.ref_seq = as_encoded_sequence_array(snps.ref_seq.ravel(), DNAArray)
    snps.alt_seq = as_encoded_sequence_array(snps.alt_seq.ravel(), DNAArray)
    kmer_indexes = get_kmer_indexes(snps.position, flank=flank)
    kmers = reference[kmer_indexes]
    forward_mask = (snps.ref_seq == "C") | (snps.ref_seq == "T")
    kmers = np.where(
        forward_mask[:, None], kmers, reverse_compliment(kmers)
    ).view(DNAArray)
    signature_encoding = MutationTypeEncoding(flank * 2 + 1)
    get_encoded_array_class(signature_encoding)
    all_hashes = signature_encoding.encode(kmers, snps)
    n_hashes = 4 ** (flank * 2) * 6
    
    if not hasattr(snps, "genotypes"):
        return count_encoded(all_hashes)
    return EncodedCounts.concatenate([count_encoded(all_hashes, weights=genotype>0)
                                      for genotype in snps.genotypes.T])

    return np.array(
        [
            np.bincount(
                all_hashes, weights=snps.genotypes[:, sample] > 0, minlength=n_hashes
            )
            for sample in range(snps.genotypes.shape[-1])
        ],
        dtype=int,
    )
