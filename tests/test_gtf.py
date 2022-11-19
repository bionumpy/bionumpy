import pytest
import bionumpy as bnp
from bionumpy import as_encoded_array
from bionumpy.gtf import get_transcript_sequences, get_exons


@pytest.fixture
def gtf_entries():
    return bnp.open("example_data/small.gtf").read()


@pytest.fixture
def reference_sequences():
    return as_encoded_array("A" * 40000, bnp.encodings.alphabet_encoding.ACGTnEncoding)


def test_get_transcript_sequences(gtf_entries, reference_sequences):
    transcript_sequences = get_transcript_sequences(gtf_entries, reference_sequences)
    print(transcript_sequences)


def test_get_exons(gtf_entries):
    assert len(get_exons(gtf_entries)) == 3
