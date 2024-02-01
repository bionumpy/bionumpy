import dataclasses
import bionumpy as bnp
import numpy as np
from npstructures.testing import assert_raggedarray_equal
from numpy.testing import assert_array_equal
from bionumpy.bnpdataclass import BNPDataClass
import bionumpy as bnp
import bionumpy.encoded_array
import bionumpy.io.vcf_buffers
from bionumpy.bnpdataclass.bnpdataclass import narrow_type
from bionumpy.io.vcf_buffers import VCFBuffer2
import pytest

from bionumpy.datatypes import VCFEntry, VCFEntryWithGenotypes
from bionumpy.encodings.vcf_encoding import PhasedGenotypeRowEncoding, GenotypeRowEncoding, PhasedHaplotypeRowEncoding
from bionumpy.util.testing import assert_bnpdataclass_equal
from tests.util import get_file_name


def test_vcf_matrix_buffer():
    f = bnp.open("example_data/variants_with_header.vcf",
                 buffer_type=bionumpy.io.vcf_buffers.PhasedVCFMatrixBuffer)

    out = bnp.open("test1.vcf", mode="w")

    for chunk in f:
        header = chunk.get_context("header")
        assert header
        out.write(chunk)

    filestart = open('test1.vcf').read(100)
    assert filestart.startswith('#'), filestart
    # check that header was written
    chunk = bnp.open("test1.vcf").read_chunk()
    assert chunk.get_context("header") != "" and chunk.get_context("header") == header


def test_vcf_matrix_buffer_stream():
    f = bnp.open("example_data/variants_with_header.vcf",
                 buffer_type=bionumpy.io.vcf_buffers.PhasedVCFMatrixBuffer)

    out = bnp.open("test1.vcf", mode="w")
    out.write(f.read_chunks())
    # check that header was written
    chunk = bnp.open("test1.vcf").read_chunk()
    assert chunk.get_context("header") != ""

def test_context_state():
    f = bnp.open("example_data/variants_with_header.vcf").read()
    assert f.get_context("header")
    f2 = bnp.open("example_data/variants.vcf").read()
    assert not f2.get_context("header")


def test_phased_genotype_encoding():
    data = bionumpy.encoded_array.as_encoded_array([
        "0|0\t0|1\t1|0\t",
        "0|0\t0|1\t1|1\t"
    ])
    encoded = PhasedGenotypeRowEncoding.encode(data)
    correct = np.array(
        [[0, 1, 2],
        [0, 1, 3]]
    )
    assert np.all(encoded == correct)


def test_phased_haplotype_encoding():
    data = bionumpy.encoded_array.as_encoded_array([
        "0|0\t0|1\t2|0\t",
        "0|0\t4|1\t1|1\t"
    ])
    encoded = PhasedHaplotypeRowEncoding.encode(data)
    correct = np.array(
        [
            [0, 0, 0, 1, 2, 0],
            [0, 0, 4, 1, 1, 1]
        ]
    )
    assert np.all(encoded == correct)


def test_genotype_encoding():
    raw_data = [
        "0|0\t0/1\t./.\t",
        "0|0\t0/0\t1|1\t",
        "0/0\t./0\t0/.\t",
        "1/1\t./0\t0/.\t",
    ]
    data = bionumpy.encoded_array.as_encoded_array(raw_data)
    encoded = bnp.EncodedArray(GenotypeRowEncoding.encode(data), bnp.encodings.BaseEncoding)
    decoded = GenotypeRowEncoding.decode(encoded)
    decoded_string = [
        ''.join((chr(c) for c in d)) + "\t" for d in decoded
    ]
    assert decoded_string == raw_data


def test_parse_unphased_vcf():
    # example_data/variants.vcf has messy unphased and missing genotypes
    filename = "example_data/variants.vcf"
    print(open(filename).read())
    f = bnp.open(filename, buffer_type=bionumpy.io.vcf_buffers.VCFMatrixBuffer)
    data = f.read()
    # assert data.get_context("header")
    decoded = GenotypeRowEncoding.decode(data.genotypes)
    print(repr(decoded))
    with bnp.open("test.tmp", "w", buffer_type=bionumpy.io.vcf_buffers.VCFMatrixBuffer) as out_file:
        out_file.write(data)
    written = [line for line in open("test.tmp").read().split("\n") if not line.startswith('#')]
    first_types = ["1|1", "1|1", "1|1", "1|1"]
    first_written = written[0].split("\t")[-4:]
    for w, t in zip(first_written, first_types):
        assert w.startswith(t), (w, t)
    # assert all(w.startswith(t) for w, t in )
    last_written = written[1].split("\t")[-4:]
    last_types = ["0/0", "1/1", "1/1", "1/1"]

    assert all(w.startswith(t) for w, t in zip(last_written, last_types))


def test_parse_phased_vcf():
    f = bnp.open("example_data/variants_phased.vcf", buffer_type=bionumpy.io.vcf_buffers.PhasedVCFMatrixBuffer)
    data = f.read()
    data = data.genotypes.raw()
    print(data)
    assert np.all(data ==
                  [
                      [3, 3, 3, 3],
                      [0, 2, 3, 3],
                      [1, 2, 1, 3]
                  ])


def test_read_info_field():
    vcf_filename = "example_data/variants_with_header.vcf"
    f = bnp.open(vcf_filename,
                 buffer_type=bionumpy.io.vcf_buffers.PhasedVCFMatrixBuffer)
    chunk = f.read_chunk()
    assert hasattr(chunk.info, 'AC')
    assert chunk.info.AC[0] == 240
    assert chunk.info.NS[0] == 2548
    assert chunk.info.NS[1] == 2548
    assert chunk.info.EX_TARGET[0] == False

@pytest.mark.skip('missing data')
def test_read_info_field2():
    vcf_filename = "example_data/info_flag.vcf"
    f = bnp.open(vcf_filename,
                 buffer_type=bionumpy.io.vcf_buffers.PhasedVCFMatrixBuffer)
    chunk = f.read_chunk()
    assert_array_equal(chunk.info.TCGA_DRIVER, [False]*6)



# @pytest.mark.xfail
def test_read_biallelic_vcf():
    file_name = "example_data/small_phased_biallelic.vcf"
    vcf = bnp.open(file_name, buffer_type=bnp.io.vcf_buffers.PhasedHaplotypeVCFMatrixBuffer)
    for chunk in vcf.read_chunks():
        print(chunk)


@pytest.mark.xfail
def test_read_info_from_vcf():
    file = "example_data/variants_with_single_individual_genotypes_and_info.vcf"
    variants = bnp.open(file).read()

    print(variants.info)
    print(variants.info.AC[0])
    print(variants.genotypes)


@pytest.mark.skip
def test_concatenate_variants():
    file = "example_data/variants_with_single_individual_genotypes_and_info.vcf"
    f = bnp.open(file)
    chunk1 = f.read_chunk(min_chunk_size=200)
    print(len(chunk1))
    chunk2 = f.read_chunk(min_chunk_size=200)
    print(len(chunk2))

    #merged = np.concatenate([chunk1, chunk2])
    #assert len(merged) == len(chunk1) + len(chunk2)
    #print(len(merged))

    other_chunk = bnp.io.vcf_buffers.VCFEntry.from_entry_tuples(
        [
            ("chr1", 1, ".", "A", "T", ".", "PASS", "AF=0.8;AC=2;AN=2")
        ]
    )

    #merged2 = np.concatenate([chunk1[0:1], other_chunk, chunk2[1:3]])

    #print(other_chunk)


# fails because comma in info field
#@pytest.mark.xfail
@pytest.fixture
def data_with_info():
    file = "example_data/vcf_symbolic_sequences.vcf"
    data = bnp.open(file).read()
    return data


def test_vcf_with_messy_info_field(data_with_info):
    for line in data_with_info:
        print(line)
        print(line.info)


def test_toiter_with_info(data_with_info):
    for entry in data_with_info.toiter():
        # assert type(entry) == bnp.io.vcf_buffers.VCFEntry.dataclass
        assert entry.info.__class__.__name__ == 'InfoDataclass'
        assert isinstance(entry.chromosome, str)
        assert isinstance(entry.position, int)
        assert isinstance(entry.info.AC, list)
        print(entry)


def test_pandas_with_info(data_with_info):
    print(dataclasses.fields(data_with_info))
    # assert issubclass(data_with_info.info, BNPDataClass), data_with_info.info

    print(
        [(f.name, f.type) for f in dataclasses.fields(data_with_info)]
    )
    df = data_with_info.topandas()
    dc = data_with_info.from_data_frame(df)
    assert_bnpdataclass_equal(data_with_info, dc)


# @pytest.mark.skip  # .genotype not implemented
def test_read_genotype_data_from_messy_vcf():
    file_name = "example_data/polaris_small.vcf"
    data = bnp.open(file_name, buffer_type=VCFBuffer2).read()
    genotype = data.genotype[0]
    assert np.all(genotype == ["0/1", "0/1", ".", "0/1"])


def test_read_genotype_with_more_data():
    file_name = get_file_name("example_data/syndip.vcf")
    data = bnp.open(file_name, buffer_type=VCFBuffer2).read()
    genotypes = data.genotype[:4]
    assert np.all(genotypes == [['1|0'], ['1|0'], ['0|1'], ['1|0']])

def test_write_genotype():
    data = narrow_type(VCFEntryWithGenotypes, 'info', str)(
        ['chr1', 'chr2'],
        [1, 2],
        ['.', 'rs123'],
        ['A', 'C'],
        ['T', 'G'],
        ['.', '.'],
        ['PASS', 'PASS'],
        ['.', '.'],
        [['0|0', '0|1'], ['0|0', '0|1']]
    )
    with bnp.open("tmp.vcf", "w", buffer_type=VCFBuffer2) as f:
        f.write(data)
    new_data = bnp.open("tmp.vcf", buffer_type=VCFBuffer2).read().get_data_object()
    assert_bnpdataclass_equal(data, new_data)

def test_read_genotype_with_no_data():
    file_name = "example_data/variants_without_genotypes.vcf"
    data = bnp.open(file_name, buffer_type=VCFBuffer2).read()
    genotypes = data.genotype[:4]
    assert genotypes.shape == (4, 0)

def test_read_empty_vcf():
    file_name = get_file_name("example_data/empty_variants.vcf")
    data = bnp.open(file_name, buffer_type=VCFBuffer2).read()
    assert len(data) == 0
    assert data.genotype.shape[0] == 0

@pytest.mark.skip   # genotype fields not implemented
def test_read_genotype_ad_field():
    file_name = "example_data/syndip.vcf"
    data = bnp.open(file_name, buffer_type=VCFBuffer2).read()
    assert_array_equal(data[0].genotype_data.AD == [1, 1])
    # AD is variable length int, so should give ragged array?
    # genotype_data or other name
    assert_raggedarray_equal(data.genotype_data.AD[0, 1] == [[1, 1], [1, 1]])


@pytest.mark.skip   # genotype fields not implemented
def test_read_genotype_ad_field():
    file_name = "example_data/syndip.vcf"
    data = bnp.open(file_name, buffer_type=VCFBuffer2).read()
    assert_array_equal(data[0].genotype_data.AD == [1, 1])
    # AD is variable length int, so should give ragged array?
    # genotype_data or other name
    assert_raggedarray_equal(data.genotype_data.AD[0, 1] == [[1, 1], [1, 1]])


def test_read_thousand_genomes_info_field():
    data = bnp.open("example_data/thousand_genomes.vcf").read()

    assert_raggedarray_equal(
        data.info.SAS_AF[0:3],
        [[0.15],
            [0.3],
            [0.16]]
    )


def test_read_hprc_multiallelic():
    data = bnp.open("example_data/hprc_multiallelic.vcf").read()

    result = data.info.AF[0:2]
    assert_raggedarray_equal(result, [
        [0.5, 0.0277778],
        [0.527778]
    ])


def test_read_write_vcf_gives_identical_file():
    file = "example_data/variants_with_single_individual_genotypes_and_info.vcf"
    variants = bnp.open(file).read()

    with bnp.open("tmp.vcf", "w") as f:
        f.write(variants)

    correct = open(file).read()
    written = open("tmp.vcf").read()

    assert written == correct


@pytest.mark.xfail
def test_read_vcf_replace_field():
    file = "example_data/variants_with_single_individual_genotypes_and_info.vcf"
    variants = bnp.open(file).read()
    variants = bnp.replace(variants, position=np.ones_like(variants.position))

    print(variants)

    with bnp.open("tmp.vcf", "w") as f:
        f.write(variants)

    data = bnp.open("tmp.vcf").read()
    assert_array_equal(data.position, variants.position)

    print(data)


#@pytest.mark.xfail
def test_parse_vcf_that_fails():
    vcf = bnp.open(get_file_name("example_data/variants_with_af.vcf")).read()
    print(vcf)


def test_ioi():
    out_filename = "tmp_ioi.vcf"
    i = bnp.open(get_file_name("example_data/thousand_genomes.vcf"),
                 buffer_type=VCFBuffer2).read()
    i = bnp.replace(i, position=i.position + 1)
    bnp.open(out_filename, "w").write(i)
    i2 = bnp.open(out_filename, buffer_type=VCFBuffer2).read()
    assert np.all(i2.genotype == i.genotype)
    assert np.all(i2.position[:10] == i.position[:10])