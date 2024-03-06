import pytest
import bionumpy as bnp
import numpy as np
from bionumpy import as_encoded_array
from bionumpy.sequence.genes import get_transcript_sequences
from bionumpy.util.testing import assert_encoded_array_equal

@pytest.fixture
def gtf_entries(data_path ):
    return bnp.open(data_path / "small.gtf").read()

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


def test_get_gene_id(gtf_entries):
    gene_id = gtf_entries.get_transcripts().gene_id[0]
    print(gene_id)
    assert gene_id == 'ENST00000619216.1'
    # assert_encoded_array_equal(gene_id, )


def test_get_exons(gtf_entries):
    assert len(gtf_entries.get_exons()) == 3


# @pytest.mark.skip('waiting')
def test_read_gff(data_path):
    annotation = bnp.open(data_path / 'small_gff.gff3').read()
    genes = annotation.get_genes()
    assert genes[0].gene_id == 'ENSG00000290825.1'
    # assert_encoded_array_equal(genes[0].gene_id, 'ENSG00000290825.1')
    # print(genes)


def test_read_sarcer_gtf(data_path):
    annotation = bnp.open(data_path / 'sacCer3.ensGene.gtf.gz').read()
    transcripts = annotation.get_transcripts()
    assert len(transcripts) > 0
