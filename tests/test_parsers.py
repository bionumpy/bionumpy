from itertools import chain
import pytest
import numpy as np
from bionumpy.file_buffers import FastQBuffer, TwoLineFastaBuffer
from bionumpy.datatypes import Interval, SNP, SequenceEntry, Variant
from bionumpy.delimited_buffers import BedBuffer, VCFBuffer, GfaSequenceBuffer
from .buffers import fastq_buffer, twoline_fasta_buffer, bed_buffer, vcf_buffer, vcf_buffer2, gfa_sequence_buffer, combos
from bionumpy.parser import chunk_lines

np.seterr(all='raise')


def chunk_from_text(text):
    return np.frombuffer(bytes(text, encoding="utf8"), dtype=np.uint8)


@pytest.mark.parametrize("buffer_name", ["bed", "vcf2", "vcf", "fastq", "fasta", "gfa_sequence"])
def test_buffer_read(buffer_name):
    buf, true_data, buf_type = combos[buffer_name]
    data = buf_type.from_raw_buffer(buf).get_data()
    assert list(data) == true_data


@pytest.mark.parametrize("buffer_name", ["fasta", "fastq"])# "bed", "vcf2", "vcf", "fastq", "fasta"])
def test_buffer_write(buffer_name):
    true_buf, data, buf_type = combos[buffer_name]
    data = data[0].stack_with_ragged(data)
    buf = buf_type.from_data(data)
    print(buf)
    print(true_buf)
    assert np.all(true_buf == buf)


def test_twoline_fasta_buffer(twoline_fasta_buffer):
    buf = TwoLineFastaBuffer.from_raw_buffer(twoline_fasta_buffer)
    seqs = buf.get_sequences()
    assert seqs.to_sequences() == ["CTTGTTGA", "CGG"]


def test_fastq_buffer(fastq_buffer):
    buf = FastQBuffer.from_raw_buffer(fastq_buffer)
    seqs = buf.get_sequences()
    assert seqs.to_sequences() == ["CTTGTTGA", "CGG"]


def test_gfa_sequence_buffer(gfa_sequence_buffer):
    buf = GfaSequenceBuffer.from_raw_buffer(gfa_sequence_buffer)
    entries = list(buf.get_sequences())
    true = [
        SequenceEntry("id1", "AACCTTGG"),
        SequenceEntry("id4", "ACTG")
    ]
    print("NAM")
    print(entries[0].name)
    assert entries == true


def test_bed_buffer(bed_buffer):
    buf = BedBuffer.from_raw_buffer(bed_buffer)
    intervals = list(buf.get_intervals())
    assert intervals == [
        Interval("chr1", 1, 3),
        Interval("chr1", 40, 60),
        Interval("chr2",  400, 600)]


def test_vcf_buffer(vcf_buffer):
    buf = VCFBuffer.from_raw_buffer(vcf_buffer)
    snps = list(buf.get_snps())
    true = [SNP("chr1",	88361,	"A",	"G"),
            SNP("chr1",	887559,	"A",	"C"),
            SNP("chr2",	8877,	"A",	"G")]
    print(true)
    print(snps)
    assert snps == true


def test_vcf_buffer2(vcf_buffer2):
    buf = VCFBuffer.from_raw_buffer(vcf_buffer2)
    variants = buf.get_variants()
    print(variants)
    true = [SNP("chr1",	88361,	"A",	"G"),
            SNP("chr1",	887559,	"A",	"CAA"),
            SNP("chr2",	8877,	"AGG",	"C")]
    assert list(variants) == true


def test_line_chunker(vcf_buffer2):
    lines = list(chain.from_iterable(chunk_lines([VCFBuffer.from_raw_buffer(vcf_buffer2).get_data()], n_lines=1)))
    true = [Variant("chr1",	88361,	"A",	"G"),
            Variant("chr1",	887559,	"A",	"CAA"),
            Variant("chr2",	8877,	"AGG",	"C")]
    print(VCFBuffer.from_raw_buffer(vcf_buffer2).get_data())
    for line, t in zip(lines, true):
        assert line == t
        
        
    

