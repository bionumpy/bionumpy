import logging
from typing import List

import numpy as np
from bionumpy.datatypes import VCFEntry
from bionumpy.genomic_data import GenomicSequence
from .sequences import simulate_sequences
from bionumpy import EncodedRaggedArray, EncodedArray, DNAEncoding
from bionumpy.encodings import ACGTnEncoding
import bionumpy as bnp
from ..bnpdataclass import bnpdataclass
from ..string_array import StringArray


@bnpdataclass
class VCFEntry:
    chromosome: str
    position: int
    id: str
    ref_seq: str
    alt_seq: str
    quality: str
    filter: str
    info: str

@bnpdataclass
class VCFEntryWithGenotypes:
    chromosome: str
    position: int
    id: str
    ref_seq: str
    alt_seq: str
    quality: str
    filter: str
    info: str
    genotypes: List[str]


def simulate_variants(genome: GenomicSequence, snp_prob=0.001, small_indel_prob=0.0001, sv_prob=0.00005,
                      ignore_variants_with_n=True, rng=np.random.default_rng(), genotype_probs={'0/0': 0.25,
                                                                                                '0/1': 0.25,
                                                                                                '1/0': 0.25,
                                                                                                '1/1': 0.25}, n_samples=0):
    chromosomes = genome.genome_context.chrom_sizes
    variant_id_offset = 0

    for chromosome, chromosome_size in chromosomes.items():
        # find lengths of variants
        lengths = np.concatenate([
            np.ones(int(chromosome_size * snp_prob)),  # SNPs
            np.random.randint(3, 50, size=int(small_indel_prob * chromosome_size)),  # small indels
            np.random.randint(50, 500, size=int(sv_prob * chromosome_size))  # SVs
        ])

        too_large = lengths > chromosome_size
        if np.any(too_large):
            logging.warning("Genome is too small to simulate some variants. Lengths are decreased")
            lengths[too_large] = chromosome_size - 3

        n_variants = len(lengths)
        logging.info("Simulating %d variants on chromosome %s", n_variants, chromosome)

        positions = rng.integers(1, chromosome_size - lengths - 1, size=n_variants)
        _, unique = np.unique(positions, return_index=True)
        positions = positions[unique]
        lengths = lengths[unique]
        n_variants = len(positions)
        is_insertion = rng.choice([True, False], n_variants)

        first_ref_base = genome[chromosome][positions]

        ref_lengths = lengths.copy()
        ref_lengths[is_insertion] = 1
        alt_lengths = lengths.copy()
        alt_lengths[~is_insertion] = 1

        # simulate alt sequences
        total_alt_length = int(sum(alt_lengths))
        alt_sequences = EncodedRaggedArray(EncodedArray(rng.integers(0, 4, size=total_alt_length), ACGTnEncoding),
                                           alt_lengths)

        ref_sequences = bnp.ragged_slice(genome[chromosome], positions, positions + ref_lengths)

        assert len(alt_sequences) == len(ref_sequences)

        # first alt sequences is always the ref base, except for SNPs
        not_snp = lengths != 1
        alt_sequences[not_snp, 0] = first_ref_base[not_snp]
        ref_sequences[:, 0] = first_ref_base

        # make sure alt sequences for SNPs do not match the ref base
        new_snp_bases = bnp.EncodedArray(
            (ref_sequences[~not_snp, 0].raw() + rng.integers(1, 4, size=sum(~not_snp))) % 4, DNAEncoding)
        new_snp_bases = bnp.change_encoding(new_snp_bases, ACGTnEncoding)
        alt_sequences[~not_snp, 0] = new_snp_bases
        assert np.all(new_snp_bases != ref_sequences[~not_snp, 0])

        variants = VCFEntry(
            chromosome=bnp.as_encoded_array([chromosome] * n_variants),
            position=positions,
            id=bnp.as_encoded_array([f"simulated{i + variant_id_offset}" for i in range(n_variants)]),
            ref_seq=bnp.change_encoding(ref_sequences, bnp.BaseEncoding),
            alt_seq=bnp.change_encoding(alt_sequences, bnp.BaseEncoding),
            quality=bnp.as_encoded_array(["."] * n_variants),
            filter=bnp.as_encoded_array(["PASS"] * n_variants),
            info=bnp.as_encoded_array(["."] * n_variants))

        if ignore_variants_with_n:
            has_n = np.any(variants.ref_seq == "N", axis=1) + np.any(variants.alt_seq == "N", axis=1)
            logging.info(f"Skipped {np.sum(has_n)} variants with N")
            variants = variants[~has_n]
        else:
            logging.info("Not skipping variants with N")

        sorting = np.argsort(variants.position)
        variants= variants[sorting]
        if n_samples == 0:
            yield variants
        else:
            n_genotypes= len(variants)*n_samples
            genotypes = np.random.choice(list(genotype_probs), size=n_genotypes, p=list(genotype_probs.values()))
            genotypes = genotypes.reshape(len(variants), n_samples)
            genotypes = StringArray(genotypes)
            yield VCFEntryWithGenotypes(
                chromosome=variants.chromosome,
                position=variants.position,
                id=variants.id,
                ref_seq=variants.ref_seq,
                alt_seq=variants.alt_seq,
                quality=variants.quality,
                filter=variants.filter,
                info=variants.info,
                genotypes=genotypes
            )

