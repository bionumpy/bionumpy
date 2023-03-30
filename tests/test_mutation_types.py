import pytest
from bionumpy.datatypes import Variant
from bionumpy.variants.mutation_signature import count_mutation_types, count_mutation_types_genomic
from bionumpy.genomic_data.genome import Genome
from bionumpy.genomic_data.genomic_sequence import GenomicSequenceDict


@pytest.fixture
def snps():
    return Variant(["chr1"]*3, [2, 3, 6], ["A", "C", "G"], ["T", "G", "A"])


@pytest.fixture
def reference():
    return "CCACCCGT"


@pytest.fixture
def reference_sequence(reference):
    return GenomicSequenceDict({'chr1': reference})


@pytest.fixture
def snp_locations(snps, reference):
    return Genome({'chr1': len(reference)}).get_locations(snps)


def test_genomic_verson(snp_locations, reference_sequence):
    counts = count_mutation_types_genomic(snp_locations, reference_sequence)
    for t in ["G[T>A]G", "A[C>G]C", "A[C>T]G"]:
        assert counts[t] == 1
    assert counts.counts.sum() == 3


def test_count_mutation_types(snps, reference):
    counts = count_mutation_types(snps, reference)
    print(counts)
    for t in ["G[T>A]G", "A[C>G]C", "A[C>T]G"]:
        assert counts[t] == 1
    assert counts.counts.sum() == len(snps)
