from bionumpy.io.regexp import match_regexp
from bionumpy.util.testing import assert_encoded_raggedarray_equal
import bionumpy as bnp
import pytest

@pytest.fixture
def search_string():
    return bnp.as_encoded_array(['''gene_id "Q0297"; transcript_id "Q0297_mRNA"; exon_number "1"; exon_id "Q0297_mRNA.1"; gene_name "Q0297";''',
                                 '''gene_id "Q0285"; transcript_id "Q0285_ncRNA"; exon_number "1"; exon_id "Q0285_ncRNA.1"; gene_name "Q0285";'''])


@pytest.fixture
def regexp():
    return r'''exon_id \"(.*?)\"'''


def test_regexp_find(search_string, regexp):
    assert_encoded_raggedarray_equal(match_regexp(search_string.ravel(), regexp), ['Q0297_mRNA.1', 'Q0285_ncRNA.1'])
