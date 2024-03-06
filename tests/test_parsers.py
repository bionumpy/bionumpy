import os
from itertools import chain
import pytest
import numpy as np
from npstructures.npdataclasses import shallow_tuple

from bionumpy.datatypes import Interval, SNP, SequenceEntry, VCFEntry
from bionumpy.io.delimited_buffers import BedBuffer, GfaSequenceBuffer, get_bufferclass_for_datatype
from bionumpy.io import VCFBuffer, FastQBuffer, TwoLineFastaBuffer
from bionumpy.io.files import bnp_open
from bionumpy.util.testing import assert_bnpdataclass_equal, assert_encoded_raggedarray_equal
from .buffers import fastq_buffer, twoline_fasta_buffer, bed_buffer, vcf_buffer, vcf_buffer2, gfa_sequence_buffer, combos, data
from bionumpy.io.parser import chunk_lines
from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.encoded_array import from_encoded_array
import bionumpy as bnp
import glob

np.seterr(all='raise')


def chunk_from_text(text):
    return np.frombuffer(bytes(text, encoding="utf8"), dtype=np.uint8)


@pytest.mark.parametrize("buffer_name", ["bed", "vcf2", "vcf", "fastq", "fasta", "gfa_sequence"])
def test_buffer_read(buffer_name):
    buf, true_data, buf_type = combos[buffer_name]
    data = buf_type.from_raw_buffer(buf).get_data()
    for line, true_line in zip(data, true_data):
        assert_bnpdataclass_equal(line, true_line)


#@pytest.mark.parametrize("buffer_name", ["fastq", "fasta", "multiline_fasta"])  # "bed", "vcf2", "vcf", "fastq", "fasta"])
@pytest.mark.parametrize("buffer_name", ["fastq"])  # "bed", "vcf2", "vcf", "fastq", "fasta"])
def test_buffer_write(buffer_name):
    true_buf, data, buf_type = combos[buffer_name]
    if buffer_name == "multiline_fasta":
        buf_type.n_characters_per_line = 6
    shallow_tuples = [shallow_tuple(d) for d in data]
    data = buf_type.dataclass.from_entry_tuples(shallow_tuples)
    # data = np.concatenate(data)
    print("DATA")
    print(repr(data))
    buf = buf_type.from_data(data)
    print("Buffer")
    print(repr(buf))
    print(buf.to_string())
    print(true_buf.to_string())
    print("BUF")
    print(repr(buf))
    print("TRUE BUF")
    print(repr(true_buf))
    print("BUF qual")
    assert np.all(true_buf == buf)


@pytest.mark.parametrize("file", ["reads.fq", "big.fq.gz"])
@pytest.mark.parametrize("chunk_size", [100, 5000000])
def test_buffered_writer_ctx_manager(file, chunk_size, tmp_path, data_path):
    file = data_path / file
    file_path = tmp_path / "tmp.fq"
    true_stream = bnp_open(data_path /'reads.fq').read_chunks()

    with bnp_open(file_path, mode='w') as f:
        f.write(true_stream)

    true_stream = bnp_open(data_path / 'reads.fq').read_chunks()
    fq_stream = bnp_open(file_path)
    for fq_item, true_item in zip(fq_stream, true_stream):
        assert_bnpdataclass_equal(fq_item, true_item)

    # os.remove(file_path)


def test_custom_read():

    @bnpdataclass
    class SampleDC:
        sequence_aa: str
        sequence: str

    for extension, delimiter in {"tsv": "\t", "csv": ","}.items():
        path = f"./tmp.{extension}"
        with open(path, 'w') as file:
            file.writelines(f"sequence{delimiter}sequence_aa\nAACCTAGGC{delimiter}ATF\nAACCTAGGC{delimiter}ATF")

        data = bnp_open(path, buffer_type=get_bufferclass_for_datatype(SampleDC, delimiter=delimiter, has_header=True)).read()
        assert [s.to_string() for s in data.sequence] == ["AACCTAGGC", "AACCTAGGC"]
        assert [s.to_string() for s in data.sequence_aa] == ["ATF", "ATF"]

        os.remove(path)


def test_raises_error_for_unsupported_types():
    with pytest.raises(RuntimeError):
        bnp_open("tmp.airr")

    with pytest.raises(RuntimeError):
        bnp_open('tmp.csv')


def test_twoline_fasta_buffer(twoline_fasta_buffer):
    buf = TwoLineFastaBuffer.from_raw_buffer(twoline_fasta_buffer)
    seqs = buf.get_data()
    assert from_encoded_array(seqs.sequence) == ["CTTGTTGA", "CGG"]


def test_fastq_buffer(fastq_buffer):
    buf = FastQBuffer.from_raw_buffer(fastq_buffer)
    seqs = buf.get_data()
    assert from_encoded_array(seqs.sequence) == ["CTTGTTGA", "CGG"]


@pytest.mark.skip("makingtrouble")
def test_read_example_data(file_name):
    if "broken" in file_name:
        # broken data should not pass tests
        return

    if file_name.endswith(".txt"):
        return

    if file_name.endswith(".sam"):
        return

    if file_name.endswith(".fa.fai"):
        return

    if file_name.endswith(".jaspar"):
        return

    if file_name.endswith(".tsv"):
        return

    if file_name.startswith("demultiplex"):
        return

    file = bnp.open(file_name)
    for chunk in file.read_chunks():
        continue

    assert True

def test_gfa_sequence_buffer(gfa_sequence_buffer):
    buf = GfaSequenceBuffer.from_raw_buffer(gfa_sequence_buffer)
    entries = list(buf.get_data())
    true = [
        SequenceEntry.single_entry("id1", "AACCTTGG"),
        SequenceEntry.single_entry("id4", "ACTG")
    ]
    for entry, t in zip(entries, true):
        assert_bnpdataclass_equal(entry, t)

@pytest.mark.skip("Replaced")
def test_vcf_buffer(vcf_buffer):
    buf = VCFBuffer.from_raw_buffer(vcf_buffer)
    snps = list(buf.get_snps())
    true = [SNP("chr1", 88361, "A", "G"),
            SNP("chr1", 887559, "A", "C"),
            SNP("chr2", 8877, "A", "G")]
    print(true)
    print(snps)
    assert snps == true


@pytest.mark.skip("Replaced")
def test_vcf_buffer2(vcf_buffer2):
    buf = VCFBuffer.from_raw_buffer(vcf_buffer2)
    variants = buf.get_variants()
    print(variants)
    true = [SNP("chr1", 88361, "A", "G"),
            SNP("chr1", 887559, "A", "CAA"),
            SNP("chr2", 8877, "AGG", "C")]
    assert list(variants) == true


def test_line_chunker(vcf_buffer2):
    lines = list(chain.from_iterable(chunk_lines([VCFBuffer.from_raw_buffer(vcf_buffer2).get_data()], n_lines=1)))
    true = data["vcf2"]
    for line, t in zip(lines, true):
        assert_bnpdataclass_equal(line, t)


def test_read_chunk_after_read_chunks_returns_empty_dataclass(data_path):
    file = bnp.open(data_path / 'reads.fq')
    chunks = list(file.read_chunks())
    new_chunk = file.read_chunk()
    assert isinstance(chunks[0],
                      type(new_chunk))


def test_read_gtf(data_path):
    file = bnp.open(data_path / 'small.gtf')
    chunk = file.read_chunk()
    assert True


def test_read_bam(data_path):
    data = bnp.open(data_path / 'alignments.bam').read()
    data2 = bnp.open(data_path / 'alignments.sam').read()
    print(data)
    print(data2)
    print(data.sequence.raw(), data.sequence.shape)

    print(data2.sequence.shape)
    assert_encoded_raggedarray_equal(data.sequence, data2.sequence)
    print(data)
    n_lines = len([line for line in open(data_path / 'alignments.sam') if not line.startswith("@")])
    assert n_lines == len(data)

