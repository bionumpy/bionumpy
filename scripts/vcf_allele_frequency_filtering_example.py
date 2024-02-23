import bionumpy as bnp


def filter_file_on_allele_count(input_file: str, output_file: str, min_ac: int = 10):
    with bnp.open(output_file, "w") as output_file:
        for chunk in bnp.open(input_file).read_chunks():
            output_file.write(chunk[chunk.info.AC.ravel() >= min_ac])


def test():
    filter_file_on_allele_count("example_data/variants_with_header.vcf", "test.vcf", min_ac=1)
    assert bnp.count_entries("test.vcf") == 53


if __name__ == "__main__":
    import sys

    filter_file_on_allele_count(sys.argv[1], sys.argv[2], int(sys.argv[3]))
