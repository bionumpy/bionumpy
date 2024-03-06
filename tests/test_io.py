import os
from pathlib import Path

import pytest
import numpy as np
from io import BytesIO
import bionumpy as bnp
from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.encodings.string_encodings import StringEncoding
from bionumpy.io.files import NumpyFileReader, NpDataclassReader, NpBufferedWriter
from .buffers import buffer_texts, combos, big_fastq_text, SequenceEntryWithQuality, VCFEntry
from bionumpy.io.matrix_dump import matrix_to_csv, parse_matrix
from bionumpy.util.testing import assert_bnpdataclass_equal, assert_encoded_array_equal, \
    assert_encoded_raggedarray_equal
from numpy.testing import assert_equal

from .util import get_file_name


@pytest.mark.parametrize("file_format", combos.keys())
@pytest.mark.parametrize("chunked", [False, True])
def test_read_write_roundtrip(file_format, chunked):
    if file_format in ("multiline_fasta", "bed12", 'wig'):
        return
    _, _, buf_type = combos[file_format]
    buffer_text = buffer_texts[file_format] * 100
    input_buffer = bytes(buffer_text, encoding="utf8")
    in_obj = BytesIO(input_buffer)
    out_buffer = bytes()
    out_obj = BytesIO(out_buffer)
    reader = NpDataclassReader(NumpyFileReader(in_obj, buffer_type=buf_type))
    writer = NpBufferedWriter(out_obj, buf_type)
    if chunked:
        for chunk in reader.read_chunks(200):
            writer.write(chunk)
    else:
        writer.write(reader.read())
    print(out_obj.getvalue())
    print(input_buffer)
    assert out_obj.getvalue() == input_buffer


@pytest.fixture
def header():
    return ["A", "AB", "ABC", "ABCD", "ABCDE"]


@pytest.fixture
def row_names():
    return ['1', '2', '3', '4']


@pytest.fixture
def integers():
    return np.arange(20).reshape(4, 5)


@pytest.fixture
def matrix_text(header, integers):
    sep = '\t'
    return sep.join(header) + "\n" + "\n".join(sep.join(str(i) for i in row) for row in integers) + "\n"


@pytest.fixture
def matrix_text2(header, integers, row_names):
    sep = '\t'
    return sep.join(['type'] + header) + "\n" + "\n".join(
        sep.join([row_name] + [str(i) for i in row]) for row_name, row in zip(row_names, integers)) + "\n"


def test_matrix_to_csv(header, integers, row_names):
    header = ["A", "AB", "ABC", "ABCD", "ABCDE"]
    integers = np.arange(20).reshape(4, 5)
    text = matrix_to_csv(integers, header=header)

    assert text.to_string() == ",".join(header) + "\n" + "\n".join(
        ",".join(str(i) for i in row) for row in integers) + "\n"


def test_read_matrix(header, integers, matrix_text):
    matrix = parse_matrix(matrix_text, rowname_type=None, field_type=int)
    assert_encoded_raggedarray_equal(matrix.col_names, header)
    assert_equal(integers, matrix.data)


def test_read_matrix_with_row_names(header, integers, matrix_text2, row_names):
    matrix = parse_matrix(matrix_text2, field_type=int)
    assert_encoded_raggedarray_equal(matrix.col_names, header)
    assert_equal(integers, matrix.data)
    assert_encoded_raggedarray_equal(matrix.row_names, row_names)


@pytest.mark.parametrize("buffer_name",
                         ["bed", "vcf2", "vcf", "fastq", "fasta", "gfa_sequence", "multiline_fasta", 'wig'])
def test_buffer_read(buffer_name):
    _, true_data, buf_type = combos[buffer_name]
    text = buffer_texts[buffer_name]
    io_obj = BytesIO(bytes(text, encoding="utf8"))
    data = NpDataclassReader(NumpyFileReader(io_obj, buf_type)).read()
    for line, true_line in zip(data, true_data):
        assert_bnpdataclass_equal(line, true_line)


@pytest.mark.parametrize("buffer_name", ["bed", "vcf2", "vcf", "fastq", "fasta", "gfa_sequence", "multiline_fasta"])
@pytest.mark.parametrize("min_chunk_size", [5000000, 50])
def test_buffer_read_chunks(buffer_name, min_chunk_size):
    _, true_data, buf_type = combos[buffer_name]
    text = buffer_texts[buffer_name]
    io_obj = BytesIO(bytes(text, encoding="utf8"))
    data = np.concatenate(list(NpDataclassReader(NumpyFileReader(io_obj, buf_type)).read_chunks(min_chunk_size)))
    for line, true_line in zip(data, true_data):
        assert_bnpdataclass_equal(line, true_line)


def test_read_big_fastq():
    io_obj = BytesIO(bytes(big_fastq_text * 20, encoding="utf8"))
    for b in NpDataclassReader(NumpyFileReader(io_obj, combos["fastq"][2])).read_chunks(min_chunk_size=1000):
        print(b)


@pytest.mark.parametrize("buffer_name", ["bed", "vcf", "fastq", "fasta"])
def test_ctx_manager_read(buffer_name, tmp_path):
    file_path = tmp_path / f"./{buffer_name}_example.{buffer_name}"
    with open(file_path, "w") as file:
        file.write(buffer_texts[buffer_name])

    with bnp.open(file_path) as file:
        file.read()


@pytest.mark.parametrize("buffer_name", ["bed", "vcf", "fastq", "fasta"])
def test_append_to_file(buffer_name):
    file_path = Path(f"./{buffer_name}_example.{buffer_name}")

    _, true_data, buf_type = combos[buffer_name]
    text = buffer_texts[buffer_name]
    io_obj = BytesIO(bytes(text, encoding="ascii"))
    data = NpDataclassReader(NumpyFileReader(io_obj, buf_type)).read()

    with bnp.open(file_path, mode="w") as file:
        file.write(data)

    with bnp.open(file_path, mode="a") as file:
        file.write(data)

    with bnp.open(file_path) as file:
        data_read = file.read()

    print(data_read)

    assert len(data_read) == len(data) * 2

    os.remove(file_path)


def test_write_dna_fastq():
    _, data, buf_type = combos["fastq"]
    entry = SequenceEntryWithQuality(["name"], ["ACGT"], ["!!!!"])
    entry.sequence = bnp.as_encoded_array(entry.sequence, bnp.DNAEncoding)
    result = buf_type.from_raw_buffer(buf_type.from_data(entry)).get_data()
    print(result)
    assert np.all(entry.sequence == result.sequence)


def test_write_empty(tmp_path):
    entry = VCFEntry([], [], [], [],
                     [], [], [], [])
    filename = tmp_path / 'tmp.vcf'
    with bnp.open(filename, 'w') as f:
        f.write(entry)


def test_read_write_bool():
    @bnpdataclass
    class BoolDC:
        param: bool

    obj = BoolDC(param=[True, False, True])
    buf_type = bnp.io.delimited_buffers.get_bufferclass_for_datatype(BoolDC, delimiter='\t', has_header=True)

    with bnp.open('tmp_bool.tsv', 'w', buf_type) as f:
        f.write(obj)

    with bnp.open('tmp_bool.tsv', 'r', buf_type) as f:
        obj2 = f.read()

    assert np.array_equal(obj.param, obj2.param)

    os.remove('tmp_bool.tsv')


@pytest.fixture
def fastq_with_carriage_return_filename(tmp_path):
    text = '''\
@test_sequence_id_here\r
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\r
+\r
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65\r
'''
    filename = tmp_path/'carriage_return.fq'
    with open(filename, 'w') as file:
        file.write(text)
    return filename


@pytest.fixture
def bed_with_carriage_return_filename(tmp_path):
    text = '''\
chr1\t1\t2\r
chr2\t3\t4
'''
    filename = tmp_path / 'carriage_return.bed'
    with open(filename, 'w') as file:
        file.write(text)
    return filename


@pytest.fixture
def fasta_with_carriage_return_filename(tmp_path):
    text = '''\
>test_sequence_id_here\r
GACTG\r
>test_sequence_id_here2\r
GACTC\r
GAG\r
'''
    filename = tmp_path/'carriage_return.fa'
    with open(filename, 'w') as file:
        file.write(text)
    return filename


def test_carriage_return_fastq(fastq_with_carriage_return_filename):
    data = bnp.open(fastq_with_carriage_return_filename).read()
    assert len(data.sequence[0]) == 60
    assert len(data.quality[0]) == 60


def test_carriage_return_bed(bed_with_carriage_return_filename):
    data = bnp.open(bed_with_carriage_return_filename).read()
    assert len(data.start) == 2
    assert len(data.stop) == 2
    assert_equal(data.stop, [2, 4])


# @pytest.mark.xfail
def test_carriage_return_fasta(fasta_with_carriage_return_filename):
    entries = bnp.open(fasta_with_carriage_return_filename).read()
    assert_encoded_raggedarray_equal(entries.sequence, ['GACTG', 'GACTCGAG'])


# @pytest.mark.xfail
def test_carriage_return_fai(fasta_with_carriage_return_filename: Path):
    # remove file if it exists
    # add .fai to the end of the file
    filename = fasta_with_carriage_return_filename
    fai_filename = filename.with_suffix(filename.suffix + '.fai')
    if os.path.exists(fai_filename):
        os.remove(fai_filename)
    fai = bnp.open_indexed(filename)
    assert_encoded_array_equal(fai['test_sequence_id_here'].raw(), 'GACTG')
    assert_encoded_array_equal(fai['test_sequence_id_here2'].raw(), 'GACTCGAG')

def test_rwr_bed_with_change(tmp_path, data_path):
    file_path = tmp_path / 'tmp_rwr.bed'
    filename = data_path  / 'alignments.bed'
    data = bnp.open(filename, buffer_type=bnp.io.Bed6Buffer).read()
    data.start = data.start + 1
    data.stop = data.stop + 1
    data.chromosome = StringEncoding(['chr1', 'chr2']).encode(data.chromosome)
    data == data[::2]
    if os.path.exists(file_path):
        os.remove(file_path)
    bnp.open(file_path, 'w', buffer_type=bnp.io.Bed6Buffer).write(data)
    text = open(file_path).read()
    assert text.startswith('chr1'), text[:10]
    print(text)
    data2 = bnp.open(file_path).read()
    assert_equal(data.start, data2.start)
    assert np.all(data.chromosome == data2.chromosome)
