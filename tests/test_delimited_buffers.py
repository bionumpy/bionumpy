import pytest

from bionumpy.encodings import GenotypeRowEncoding
from bionumpy.io.delimited_buffers import BedBuffer
from bionumpy.datatypes import Interval
import numpy as np
import bionumpy as bnp
from npstructures.testing import assert_raggedarray_equal
import npstructures as nps

intervals = Interval(["chr1", "chr1"], [2, 10], [100, 20])



def test_get_column_range():
    data = """\
a b c d... e f
1 2 3... 4 5 6
1 asefasefasefsae 3... 4 5 6
1 asefasefasefsae 3... 892398 5 6
""".replace(" ", "\t")

    buffer = bnp.io.delimited_buffers.DelimitedBuffer.from_raw_buffer(
        np.array([ord(c) for c in data], dtype=np.uint8),
    )

    result = buffer.get_column_range_as_text(2, 5)
    correct = bnp.as_encoded_array([
        [c for c in "\t".join(row.split("\t")[2:5])]
         for row in data.split("\n")
    ][:-1])
    assert_raggedarray_equal(result, correct)


def test_from_data():
    buf = BedBuffer.from_data(intervals)
    assert buf.to_string() == """\
chr1\t2\t100
chr1\t10\t20
"""


def test_vcf_matrix_buffer():
    f = bnp.open("example_data/variants_with_header.vcf",
                 buffer_type=bnp.io.delimited_buffers.VCFMatrixBuffer)


    out = bnp.open("test1.vcf", mode="w")

    for chunk in f:
        #print(chunk)
        #print(chunk.genotypes)
        #string_genotypes = bnp.encodings.PhasedGenotypeEncoding.decode(chunk.genotypes)
        #print(string_genotypes)

        #print(chunk)
        out.write(chunk)



def test_genotype_encoding():
    data = bnp.as_encoded_array([
        "0|0\t0|1\t1|0\t",
        "0|0\t0|1\t1|1\t"
    ])
    encoded = GenotypeRowEncoding.encode(data)
    correct = np.array(
        [[0, 2, 1],
        [0, 2, 3]]
    )
    assert np.all(encoded == correct)

    #decoded = GenotypeRowEncoding.decode(encoded)
    #print(repr(decoded))
    #print(repr(data))
    #assert assert_raggedarray_equal(decoded, data)
    #assert False
