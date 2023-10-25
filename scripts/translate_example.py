import bionumpy as bnp
from bionumpy.sequence import translate_dna_to_protein
import time


def main(input_file, output_file):
    chunks = bnp.open(input_file).read_chunks()
    translated_chunks = (
        bnp.replace(chunk, sequence=translate_dna_to_protein(chunk.sequence))
        for chunk in chunks)
    with bnp.open(output_file, "w") as f:
        for chunk in translated_chunks:
            f.write(chunk)


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])