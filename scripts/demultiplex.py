import bionumpy as bnp
from pathlib import Path
"""
DEMULTIPLEX-ing
===============

This file is an working example of a demultiplex task via bionumpy.

File required:

Barcodes:
    'demultiple_R1_a.txt'
    'demultiple_R1_b.txt'
    'demultiple_R2.txt'

Reads:
    'example_1_R1_001.fastq'
    'example_1_R2_001.fastq'

Execute:
    python demultiplex.py

Output:
    List of files:
    'out_p_X_a.fastq' 
    'out_p_X_b.fastq' 
    
    where X in in the range of the number of given barcodes.

"""

directory = Path("example_data", "demultiplex")

file_sample_r1 = 'small_1_R1_001.fastq'
file_sample_r2 = 'small_1_R2_001.fastq'

filename_out = 'out'


def get_barcodes(barcode_filenames):
    barcodes = {name: [line.strip() for line in open(directory / file_name)]
                for file_name, name in zip(barcode_filenames, ['a1', 'a2', 'b'])}
    return barcodes


file_barcode_r1_a = 'demultiple_R1_a.txt'
file_barcode_r1_b = 'demultiple_R1_b.txt'
file_barcode_r2 = 'demultiple_R2.txt'

barcodes = get_barcodes([file_barcode_r1_a, file_barcode_r1_b, file_barcode_r2])


def handle_barcodes(i_file, i, j, bar_1, bar_2):
    chunks_1 = bnp.open(directory / file_sample_r1).read_chunks()
    chunks_2 = bnp.open(directory / file_sample_r2).read_chunks()

    with bnp.open(directory / f"{filename_out}__{i_file+1}__{i}_{j}__R1.fastq", 'w') as or1, \
         bnp.open(directory / f"{filename_out}__{i_file+1}__{i}_{j}__R2.fastq", 'w') as or2:
        for chunk1, chunk2 in zip(chunks_1, chunks_2):
            matches = [bnp.sequence.match_string(chunk1.sequence, bar_1),
                       bnp.sequence.match_string(chunk2.sequence, bar_2)]
            seq_matches = matches[0].any(axis=-1) & matches[1].any(axis=-1)
            or1.write(chunk1[seq_matches])
            or2.write(chunk2[seq_matches])


for i, j in ('a1', 'b'), ('a2', 'b'):
    for i_file, (bar_1, bar_2) in enumerate(zip(barcodes[i], barcodes[j])):
        handle_barcodes(i_file, i, j, bar_1, bar_2)
