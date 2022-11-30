import pytest
import bionumpy as bnp
import numpy as np
from bionumpy import as_encoded_array
from bionumpy.sequence.genes import get_transcript_sequences


@pytest.fixture
def gtf_entries():
    return bnp.open("example_data/small.gtf").read()


@pytest.fixture
def reference_sequences():
    return as_encoded_array("A" * 40000, bnp.encodings.BaseEncoding)


def test_get_transcript_sequences(gtf_entries, reference_sequences):
    transcript_sequences = get_transcript_sequences(gtf_entries, reference_sequences)
    print("TRUE")
    #true = bnp.as_encoded_array(['T'*(gtf_entries.stop[0]-gtf_entries.start[0]), 'A'*588], bnp.encodings.base_encoding.ASCIIEncoding)
    true = ['T'*(gtf_entries.stop[0]-gtf_entries.start[0]), 'A'*588]
    print(true)
    assert np.all(transcript_sequences.sequence==true)



def test_get_exons(gtf_entries):
    assert len(gtf_entries.get_exons()) == 3
