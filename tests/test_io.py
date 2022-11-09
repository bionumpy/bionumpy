import os
from pathlib import Path

import pytest
import numpy as np
from io import BytesIO
import bionumpy as bnp
from bionumpy import BedBuffer
from bionumpy.io.files import NumpyFileReader, NpDataclassReader
from .buffers import buffer_texts, combos, big_fastq_text
from bionumpy.io.matrix_dump import matrix_to_csv
from npstructures.testing import assert_npdataclass_equal


def test_matrix_to_csv():
    header = ["A", "AB", "ABC", "ABCD", "ABCDE"]
    integers = np.arange(20).reshape(4, 5)
    text = matrix_to_csv(integers, header=header)

    assert text.to_string() == ",".join(header) + "\n" + "\n".join(",".join(str(i) for i in row) for row in integers)+"\n"


@pytest.mark.parametrize("buffer_name", ["bed", "vcf2", "vcf", "fastq", "fasta", "gfa_sequence", "multiline_fasta"])
def test_buffer_read(buffer_name):
    _, true_data, buf_type = combos[buffer_name]
    text = buffer_texts[buffer_name]
    io_obj = BytesIO(bytes(text, encoding="utf8"))
    data = NpDataclassReader(NumpyFileReader(io_obj, buf_type)).read()
    for line, true_line in zip(data, true_data):
        print("#####", line, true_line)
        assert_npdataclass_equal(line, true_line)


@pytest.mark.parametrize("buffer_name", ["bed", "vcf2", "vcf", "fastq", "fasta", "gfa_sequence", "multiline_fasta"])
@pytest.mark.parametrize("min_chunk_size", [5000000, 50])
def test_buffer_read_chunks(buffer_name, min_chunk_size):
    _, true_data, buf_type = combos[buffer_name]
    text = buffer_texts[buffer_name]
    io_obj = BytesIO(bytes(text, encoding="utf8"))
    data = np.concatenate(list(NpDataclassReader(NumpyFileReader(io_obj, buf_type)).read_chunks(min_chunk_size)))
    for line, true_line in zip(data, true_data):
        assert_npdataclass_equal(line, true_line)


def test_read_big_fastq():
    io_obj = BytesIO(bytes(big_fastq_text*20, encoding="utf8"))
    for b in NpDataclassReader(NumpyFileReader(io_obj, combos["fastq"][2])).read_chunks(min_chunk_size=1000):
        print(b)


@pytest.mark.parametrize("buffer_name", ["bed", "vcf", "fastq", "fasta"])
def test_ctx_manager_read(buffer_name):
    file_path = Path(f"./{buffer_name}_example.{buffer_name}")

    with open(file_path, "w") as file:
        file.write(buffer_texts[buffer_name])

    with bnp.open(file_path) as file:
        file.read()

    os.remove(file_path)


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
