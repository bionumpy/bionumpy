# Various tests for reading and parsing real data files
import pytest
from npstructures.testing import assert_raggedarray_equal

import bionumpy as bnp
import numpy as np


def test_read_polaris_vcf():
    data = bnp.open("example_data/polaris.vcf")

    for chunk in data:
        print(chunk)
        print(chunk.info.NUM_ALT_KIDS)


def test_read_syndip_vcf():
    data = bnp.open("example_data/syndip.vcf").read()
    print(data.info)


def test_read_vcf_info_field_with_missing_header():
    data = bnp.open("example_data/vcf_with_broken_header.vcf").read()
