import pytest
import bionumpy as bnp
import numpy as np
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
    assert np.all(transcript_sequences.sequence==['T'*(gtf_entries.stop[0]-gtf_entries.start[0]), 'A'*588])



def test_get_exons(gtf_entries):
    assert len(get_exons(gtf_entries)) == 3
