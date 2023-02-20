import pytest
from bionumpy.genomic_data import GenomicIntervals, GenomicSequence
from bionumpy.encoded_array import as_encoded_array
from bionumpy.util.testing import assert_encoded_raggedarray_equal
from bionumpy.encodings import ACGTnEncoding


@pytest.fixture
def sequence_dict():
    return {'chr1': 'acctgg',
            'chr2': 'aggtg',
            'chr3': 'aaaa'}


@pytest.fixture
def genomic_sequence(sequence_dict):
    return GenomicSequence.from_dict(sequence_dict)


@pytest.fixture
def chrom_sizes(sequence_dict):
    return {name: len(seq) for name, seq in sequence_dict.items()}


@pytest.fixture
def genomic_intervals(chrom_sizes):
    return GenomicIntervals.from_fields(chrom_sizes,
                                        ['chr1', 'chr2', 'chr2'],
                                        [0, 2, 3],
                                        [2, 4, 4],
                                        ['-', '+', '-'])


@pytest.fixture
def interval_sequences():
    return as_encoded_array(['gt', 'gt', 'a'], ACGTnEncoding)


def test_extract_intervals(genomic_sequence: GenomicSequence, genomic_intervals, interval_sequences):
    assert_encoded_raggedarray_equal(
        genomic_sequence[genomic_intervals],
        interval_sequences)
