import pytest
import bionumpy as bnp
from bionumpy import as_encoded_array
from bionumpy.util.testing import assert_encoded_array_equal
from bionumpy.sequence import get_strand_specific_sequences, get_reverse_complement


@pytest.fixture
def dna_sequence():
    return as_encoded_array('ACGTACGTACGT', bnp.DNAEncoding)


@pytest.fixture
def ascii_sequence():
    return as_encoded_array('ACGTG')


def test_reverse_ascii_complement(ascii_sequence):
    reverse_complement = get_reverse_complement(ascii_sequence)
    assert_encoded_array_equal(reverse_complement, "CACGT")


@pytest.fixture
def stranded_intervals():
    return bnp.datatypes.Bed6(["chr1", "chr1"], [1, 4], [3, 7], [".", "."], [".", "."], ["+", "-"])


def test_strand_specific_sequences(dna_sequence, stranded_intervals):
    truth = as_encoded_array(['CG', 'CGT'], bnp.DNAEncoding)
    result = get_strand_specific_sequences(dna_sequence, stranded_intervals)
    bnp.util.testing.assert_encoded_raggedarray_equal(truth, result)
