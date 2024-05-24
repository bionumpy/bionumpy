import pytest

from bionumpy import MultiLineFastaBuffer
import bionumpy as bnp
from bionumpy.alignments.msa import MultipleSequenceAlignment


@pytest.fixture
def mfa_obj(data_path):
    return bnp.open(data_path / "example.mfa", buffer_type=MultiLineFastaBuffer).read()


@pytest.fixture()
def msa(mfa_obj):
    return MultipleSequenceAlignment.from_sequence_entries(mfa_obj)


def test_from_sequence_entries(mfa_obj):
    sequences = mfa_obj
    msa = MultipleSequenceAlignment.from_sequence_entries(sequences)
    msa.matrix.shape == (3, 5)

@pytest.mark.skip
def test_mask(msa):
    assert msa.mask().shape == msa.matrix.shape
    assert msa.mask().dtype == bool
    assert msa.mask().sum() == 20
