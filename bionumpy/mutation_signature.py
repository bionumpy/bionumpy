import numpy as np
from .encodings import ACTGEncoding
from .sequences import as_sequence_array
from .dna import complement
from .chromosome_map import ChromosomeMap
import logging

logger = logging.getLogger(__name__)


def get_kmer_indexes(position, flank=2):
    return position[..., np.newaxis] + np.arange(-flank, flank + 1)


class SNPEncoding:
    lookup = np.zeros((256, 256), dtype=np.uint8)
    lookup[ord("C")][[ord(c) for c in "AGT"]] = np.arange(3)
    lookup[ord("G")][[ord(c) for c in "TCA"]] = np.arange(3)
    lookup[ord("T")][[ord(c) for c in "ACG"]] = 3 + np.arange(3)
    lookup[ord("A")][[ord(c) for c in "TGC"]] = 3 + np.arange(3)

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
    def __init__(self, k, encoding=ACTGEncoding):
        self.k = k
        self.h = 4 ** np.arange(k)
        self.h[k // 2 + 1 :] = self.h[k // 2 : -1]
        self.h[k // 2] = 0
        self._encoding = encoding

    def encode(self, kmer, snp):
        kmer = as_sequence_array(kmer, self._encoding)
        assert kmer.shape[-1] == self.k, (kmer.shape, self.k)
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
    assert np.all(reference[snps.position] == snps.ref_seq)
    kmer_indexes = get_kmer_indexes(snps.position, flank=flank)
    kmers = reference[kmer_indexes]
    forward_mask = (snps.ref_seq == ord("C")) | (snps.ref_seq == ord("T"))
    kmers = np.where(
        forward_mask[:, None], kmers, complement(kmers[:, ::-1])
    )
    signature_encoding = MutationSignatureEncoding(flank * 2 + 1)
    all_hashes = signature_encoding.encode(kmers, snps)
    n_hashes = 4 ** (flank * 2) * 6

    # counter = EncodedCounter(signature_encoding)
    # if not hash(snps, "genotypes"):
    #     counter.count(all_hashes)
    # else:
    #     counter.count(all_hashes, weights=snps.genotypes.T > 0)

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
