import pytest
import numpy as np
from bionumpy.kmers import TwoBitHash
from bionumpy.file_buffers import FastQBuffer, TwoLineFastaBuffer
from bionumpy.bed_parser import Interval, SNP
from bionumpy.delimited_buffers import BedBuffer, VCFBuffer
np.seterr(all='raise')
from .buffers import fastq_buffer, twoline_fasta_buffer, bed_buffer, vcf_buffer, vcf_buffer2, combos

def chunk_from_text(text):
    return np.frombuffer(bytes(text, encoding="utf8"), dtype=np.uint8)


@pytest.mark.parametrize("buffer_name", ["bed", "vcf2", "vcf", "fastq", "fasta"])
def test_buffer_read(buffer_name):
    buf, true_data, buf_type = combos[buffer_name]
    data = buf_type.from_raw_buffer(buf).get_data()
    assert list(data) == true_data

@pytest.mark.parametrize("buffer_name", ["fasta"])# "bed", "vcf2", "vcf", "fastq", "fasta"])
def test_buffer_write(buffer_name):
    true_buf, data, buf_type = combos[buffer_name]
    data = data[0].stack_with_ragged(data)
    buf = buf_type.from_data(data)
    assert np.all(true_buf == buf)

def test_twoline_fasta_buffer(twoline_fasta_buffer):
    buf = TwoLineFastaBuffer.from_raw_buffer(twoline_fasta_buffer)
    seqs = buf.get_sequences()
    assert seqs.to_sequences() == ["CTTGTTGA", "CGG"]

def test_fastq_buffer(fastq_buffer):
    buf = FastQBuffer.from_raw_buffer(fastq_buffer)
    seqs = buf.get_sequences()
    assert seqs.to_sequences() == ["CTTGTTGA", "CGG"]

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
