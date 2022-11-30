import numpy as np

import bionumpy as bnp
import bionumpy.encoded_array
from bionumpy.encodings.vcf_encoding import PhasedGenotypeRowEncoding, GenotypeRowEncoding


def test_vcf_matrix_buffer():
    f = bnp.open("example_data/variants_with_header.vcf",
                 buffer_type=bnp.io.delimited_buffers.PhasedVCFMatrixBuffer)


    out = bnp.open("test1.vcf", mode="w")

    for chunk in f:
        #print(chunk)
        #print(chunk.genotypes)
        #string_genotypes = bnp.encodings.PhasedGenotypeEncoding.decode(chunk.genotypes)
        #print(string_genotypes)

        #print(chunk)
        out.write(chunk)


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
    f = bnp.open("example_data/variants.vcf", buffer_type=bnp.VCFMatrixBuffer)
    data = f.read()
    decoded = GenotypeRowEncoding.decode(data.genotypes)
    print(repr(decoded))
    out_file = bnp.open("test.tmp", "w", buffer_type=bnp.VCFMatrixBuffer)
    out_file.write(data)

    written = open("test.tmp").read().split("\n")
    assert written[0].split("\t")[-4:] == ["1|1", "1|1", "1|1", "1|1"]
    assert written[1].split("\t")[-4:] == ["0/0", "1/1", "1/1", "1/1"]


def test_parse_phased_vcf():
    f = bnp.open("example_data/variants_phased.vcf", buffer_type=bnp.PhasedVCFMatrixBuffer)
    data = f.read()
    data = data.genotypes.raw()
    print(data)
    assert np.all(data ==
                  [
                      [3, 3, 3, 3],
                      [0, 2, 3, 3],
                      [1, 2, 1, 3]
                  ])
