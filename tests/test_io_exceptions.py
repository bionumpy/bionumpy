from bionumpy.io.file_buffers import FastQBuffer, TwoLineFastaBuffer
from bionumpy.io.delimited_buffers import BedBuffer
from bionumpy.io.exceptions import FormatException
import bionumpy as bnp
import pytest

malformed_fastqs = ["""\
@header
actg
-
!!!!
""",
                    """\
header
actg
+
!!!!
""",
                    """\
@header
actg
+
@header
actg
+
@header
actg
+
"""]

malformed_two_line_fastas = ["""\
>header
acggtt
acggtt
>header
acgtt
"""]


malformed_bed_files = ["""\
chr1\t10\t20
chr2\t10\ttwenty
""",
                       """\
chr1\t10\t20
chr2\t10\t20\t30
"""]


@pytest.mark.parametrize("text", malformed_fastqs)
def test_fastq_raises_format_exception(text):
    buf_type = FastQBuffer
    with pytest.raises(FormatException):
        buf = buf_type.from_raw_buffer(bnp.as_encoded_array(text))
        buf.get_data()


@pytest.mark.parametrize("text", malformed_two_line_fastas)
def test_fasta_raises_format_exception(text):
    buf_type = TwoLineFastaBuffer
    with pytest.raises(FormatException):
        buf = buf_type.from_raw_buffer(bnp.as_encoded_array(text))
        buf.get_data()


@pytest.mark.parametrize("text", malformed_bed_files)
def test_bed_raises_format_exception(text):
    buf_type = BedBuffer
    with pytest.raises(FormatException):
        buf = buf_type.from_raw_buffer(bnp.as_encoded_array(text))
        buf.get_data()
