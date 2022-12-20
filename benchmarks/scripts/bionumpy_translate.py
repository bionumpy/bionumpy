import bionumpy as bnp
from bionumpy.sequence import translate_dna_to_protein
import time

t0 = time.perf_counter()
input_stream = bnp.open(snakemake.input[0]).read_chunks()
output_stream = bnp.open(snakemake.output[0], "w")
output_stream.write(translate_dna_to_protein(input_stream))
output_stream.close()

print("Translate took", (time.perf_counter()-t0))