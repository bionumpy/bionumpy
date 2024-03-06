import pytest
import numpy as np
from bionumpy.datatypes import Variant
from bionumpy.variants import count_mutation_types, count_mutation_types_genomic, MutationTypeEncoding
from bionumpy.variants.mutation_signature import MutationTypeEncoding
from bionumpy.genomic_data.genome import Genome
from bionumpy.genomic_data.genomic_sequence import GenomicSequenceDict
import bionumpy as bnp

@pytest.fixture
def snps():
    return Variant(["chr1"]*3, [2, 3, 6], ["A", "C", "G"], ["T", "G", "A"])


@pytest.fixture
def reference():
    return "CCACCCGT"


@pytest.fixture
def mutation_encoding():
    return MutationTypeEncoding(1)


@pytest.fixture
def reference_sequence(reference):
    return GenomicSequenceDict({'chr1': reference})


@pytest.fixture
def snp_locations(snps, reference):
    return Genome({'chr1': len(reference)}).get_locations(snps)


def test_cosmic_read(data_path):
    matrix = bnp.io.read_matrix(data_path / 'COSMIC_v3.3.1_SBS_GRCh38.txt')
    encoded = bnp.as_encoded_array(matrix.row_names.to_numpy_array(),
                                   MutationTypeEncoding(1))
    matrix.data[np.argsort(encoded)]

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


@pytest.mark.parametrize('s', ['A[C>A]A', 'C[T>G]A'])
def test_mutation_encoding(s, mutation_encoding):
    encoded = mutation_encoding.encode(bnp.as_encoded_array(s))
    decoded = encoded.encoding.decode(encoded.raw())
    assert decoded == s, (decoded, s)
