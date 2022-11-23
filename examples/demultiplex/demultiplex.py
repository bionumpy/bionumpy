import numpy as np
import bionumpy as bnp
from collections import defaultdict

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

file_barcode_r1_a = 'demultiple_R1_a.txt'
file_barcode_r1_b = 'demultiple_R1_b.txt'
file_barcode_r2 = 'demultiple_R2.txt'

file_sample_r1 = 'example_1_R1_001.fastq'
file_sample_r2 = 'example_1_R2_001.fastq'

filename_out = 'out'


class MatchSequence(bnp.rollable.RollableFunction):
    def __init__(self, matching_sequence):
        self._matching_sequence = bnp.sequence.as_sequence_array(matching_sequence)
        self.window_size = len(matching_sequence)

    def __call__(self, sequence):
        return np.all(sequence == self._matching_sequence, axis=-1)


barcoders = defaultdict(dict)
for file, name in zip([file_barcode_r1_a, file_barcode_r1_b, file_barcode_r2], ['a1', 'a2', 'b']):
    with open(file) as i_f:
        for lines in i_f:
            barcode = lines.strip()
            barcoders[name][barcode] = MatchSequence(barcode)

for i, j in ('a1', 'b'), ('a2', 'b'):
    for i_file, (bar_1, bar_2) in enumerate(zip(barcoders[i].values(), barcoders[j].values())):
        with bnp.open(f"{filename_out}__{i_file+1}__{i}_{j}__R1.fastq", 'w') as or1, \
                bnp.open(f"{filename_out}__{i_file+1}__{i}_{j}__R2.fastq", 'w') as or2:
            with bnp.open(file_sample_r1) as source_R1, bnp.open(file_sample_r2) as source_R2:
                for chunk1, chunk2 in zip(source_R1, source_R2):
                    matches = [bar_1.rolling_window(chunk1.sequence),
                               bar_2.rolling_window(chunk2.sequence)]
                    seq_matches = matches[0].any(axis=-1) * matches[1].any(axis=-1)
                    or1.write(chunk1[seq_matches])
                    or2.write(chunk2[seq_matches])

