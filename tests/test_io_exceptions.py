from bionumpy.io.file_buffers import FastQBuffer, TwoLineFastaBuffer
from bionumpy.io.delimited_buffers import BedBuffer
from bionumpy.io.exceptions import FormatException
import bionumpy as bnp
import pytest

malformed_fastqs = [("""\
@header
actg
-
!!!!
""", 2),
                    ("""\
header
actg
+
!!!!
""", 0),
                    ("""\
@header
actg
+
@header
actg
+
@header
actg
+
""", 4)]

malformed_two_line_fastas = [("""\
>header
acggtt
acggtt
>header
acgtt
""", 2)]


malformed_bed_files = [("""\
chr1\t10\t20
chr2\t10\ttwenty
""", 1),
                       ("""\
chr1\t10\t20
chr2\t10\t20\t30
""", 1)]


@pytest.mark.parametrize("data", malformed_fastqs)
def test_fastq_raises_format_exception(data):
    text, error_line = data
    buf_type = FastQBuffer
    with pytest.raises(FormatException) as e:
        buf = buf_type.from_raw_buffer(bnp.as_encoded_array(text))
        buf.get_data()
    assert e.value.line_number == error_line


@pytest.mark.parametrize("data", malformed_two_line_fastas)
def test_fasta_raises_format_exception(data):
    text, error_line = data
    buf_type = TwoLineFastaBuffer
    with pytest.raises(FormatException) as e:
        buf = buf_type.from_raw_buffer(bnp.as_encoded_array(text))
        buf.get_data()
    assert e.value.line_number == error_line


@pytest.mark.parametrize("data", malformed_bed_files)
def test_bed_raises_format_exception(data):
    text, error_line = data
    buf_type = BedBuffer
    with pytest.raises(FormatException) as e:
        buf = buf_type.from_raw_buffer(bnp.as_encoded_array(text))
        buf.get_data()
    # assert e.value.line_number == error_line
