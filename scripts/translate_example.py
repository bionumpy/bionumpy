import bionumpy as bnp
from bionumpy.sequence import translate_dna_to_protein


def translate_dna(input_file: str, output_file: str):
    with bnp.open(output_file, "w") as f:
        for chunk in bnp.open(input_file).read_chunks():
            translated = bnp.replace(chunk, sequence=translate_dna_to_protein(chunk.sequence))
            f.write(translated)


def test_translate_dna():
    input_file = "example_data/dna_translatable.fa"
    output_file = "tmp.fa"
    translate_dna(input_file, output_file)
    assert bnp.count_entries(output_file) == bnp.count_entries(input_file)


if __name__ == "__main__":
    import sys
    translate_dna(*sys.argv[1:])
