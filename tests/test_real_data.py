# Various tests for reading and parsing real data files
import pytest
from npstructures.testing import assert_raggedarray_equal

import bionumpy as bnp
import numpy as np

from bionumpy import EncodedRaggedArray


def test_read_polaris_vcf(data_path):
    data = bnp.open(data_path / "polaris.vcf")

    for chunk in data:
        print(chunk)
        print(chunk.info.NUM_ALT_KIDS)


def test_read_syndip_vcf(data_path):
    data = bnp.open(data_path / "syndip.vcf").read()
    print(data.info)


def test_read_vcf_info_field_with_missing_header(data_path):
    data = bnp.open(data_path / "vcf_with_broken_header.vcf").read()
    assert isinstance(data.info, EncodedRaggedArray) and data.info.encoding == bnp.BaseEncoding, \
        "Should parse as string when info tags missing"
