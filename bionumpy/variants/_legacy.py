import numpy as np
from warnings import warn
from ..genomic_data import GenomicSequence, GenomicLocation
from ..encodings import DNAEncoding
from ..datatypes import Variant
from ..encoded_array import EncodedArray
from ..encoded_array import as_encoded_array
from ..sequence import get_reverse_complement, count_encoded
from ..streams import streamable
from ..sequence.lookup import Lookup
from .mutation_signature import MutationTypeEncoding, encode_snps


def get_kmer_indexes(position, flank=1):
    return np.add.outer(position, np.arange(-flank, flank + 1))

import logging

logger = logging.getLogger(__name__)


@streamable(reduction=sum)
def count_mutation_types(variants, reference, flank=1):
    warn('count_mutation_types is depracated, use count_mutation_types_genomic instead', DeprecationWarning)
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
    hashes = encode_snps(kmers, snps.alt_seq)
    #hashes = MutationTypeEncoding(flank).encode(kmers, snps)
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
