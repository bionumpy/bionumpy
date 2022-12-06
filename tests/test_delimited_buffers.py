import bionumpy.encoded_array
from bionumpy.io.delimited_buffers import BedBuffer
from bionumpy.datatypes import Interval
import numpy as np
import bionumpy as bnp
from npstructures.testing import assert_raggedarray_equal

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
    correct = bionumpy.encoded_array.as_encoded_array([
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



