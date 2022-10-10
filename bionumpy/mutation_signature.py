import numpy as np
from .encodings import ACTGEncoding
from .sequences import as_encoded_sequence_array, as_sequence_array
from .dna import complement
from .chromosome_map import ChromosomeMap
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

    text = np.array([f"C>{c}" for c in "AGT"] + [f"T->{c}" for c in "ACG"])

    @classmethod
    def to_string(cls, encoded):
        return cls.text[encoded]

    @classmethod
    def encode(cls, snp):
        return cls.lookup[snp.ref_seq, snp.alt_seq]

    @classmethod
    def decode(cls, encoded):
        pass


class MutationSignatureEncoding:
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
        return kmer_hashes + 4 ** (self.k - 1) * snp_hashes

    def decode(self, encoded):
        snp = SNPEncoding.decode(encoded >> (2 * (self.k - 1)))
        chars = (encoded >> (2 * np.arange(self.k - 1))) & 3
        kmer = "".join(chr(b) for b in ACTGEncoding.decode(chars))
        return kmer[: self.k // 2] + "[" + snp + "]" + kmer[self.k // 2 :]

    def to_string(self, encoded):
        snp = SNPEncoding.to_string(encoded >> (2 * (self.k - 1)))
        chars = (encoded >> (2 * np.arange(self.k - 1))) & 3
        kmer = "".join(chr(b) for b in ACTGEncoding.decode(chars))
        return kmer[: self.k // 2] + "[" + snp + "]" + kmer[self.k // 2 :]


@ChromosomeMap(reduction=sum)
def get_kmers(snps, reference, flank):
    reference = as_encoded_sequence_array(reference, DNAArray)
    snps.ref_seq = as_encoded_sequence_array(snps.ref_seq, DNAArray)
    snps.alt_seq = as_encoded_sequence_array(snps.alt_seq, DNAArray)
    kmer_indexes = get_kmer_indexes(snps.position, flank=flank)
    kmers = reference[kmer_indexes]
    forward_mask = (snps.ref_seq == "C") | (snps.ref_seq == "T")
    kmers = np.where(
        forward_mask[:, None], kmers, complement(kmers[:, ::-1])
    ).view(DNAArray)
    signature_encoding = MutationSignatureEncoding(flank * 2 + 1)
    all_hashes = signature_encoding.encode(kmers, snps)
    n_hashes = 4 ** (flank * 2) * 6

    if not hasattr(snps, "genotypes"):
        counts = np.bincount(all_hashes, minlength=n_hashes)
        return counts
    return np.array(
        [
            np.bincount(
                all_hashes, weights=snps.genotypes[:, sample] > 0, minlength=n_hashes
            )
            for sample in range(snps.genotypes.shape[-1])
        ],
        dtype=int,
    )
