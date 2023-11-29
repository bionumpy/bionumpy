from io import BytesIO

from bionumpy.io import FastQBuffer, TwoLineFastaBuffer
from bionumpy.io.parser import NumpyFileReader
from bionumpy.io.delimited_buffers import BedBuffer
from bionumpy.io.exceptions import FormatException, ParsingException
from bionumpy.io.files import NpDataclassReader
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
""", 1)]
#                       ("""\
# chr1\t10\t20
#chr2\t10\t20\t30
# """, 1)]


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
    print(text)
    with pytest.raises(FormatException) as e:
        buf = buf_type.from_raw_buffer(bnp.as_encoded_array(text))
        buf.get_data()
    assert e.value.line_number == error_line


@pytest.mark.parametrize("data", malformed_fastqs)
def test_npdataclass_raises_format_exception(data):
    valid_fastq = """\
@header
acgtt
+
!!!!!
"""
    malformed, line_number = data
    fobj = BytesIO(bytes(valid_fastq*100+malformed, encoding="ascii"))
    npfilereader = NumpyFileReader(fobj, buffer_type=FastQBuffer)
    reader = NpDataclassReader(npfilereader)
    with pytest.raises(FormatException) as e:
        for _ in reader.read_chunks(200):
            pass
    assert e.value.line_number == 4*100+line_number


@pytest.mark.parametrize("data", malformed_bed_files)
def test_npdataclass_raises_format_exception_bed(data):
    valid_bed = """\
chr1\t10\t20
chr2\t20\t30
chr1\t10\t20
chr2\t20\t30
"""
    malformed, line_number = data
    fobj = BytesIO(bytes(valid_bed*100+malformed, encoding="ascii"))
    npfilereader = NumpyFileReader(fobj, buffer_type=BedBuffer)
    reader = NpDataclassReader(npfilereader)
    with pytest.raises(FormatException) as e:
        for chunk in reader.read_chunks(200):
            chunk.stop
    assert e.value.line_number == 4*100+line_number

