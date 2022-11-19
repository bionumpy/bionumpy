import bionumpy as bnp

fastq_entries_stream = bnp.open("example_data/reads.fq")
for fastq_entries in fastq_entries_stream:
    print(fastq_entries)

variants_stream = bnp.open("example_data/variants.vcf")
for chromosome, variants in variants_stream:
    print(variants)


import numpy as np
sequences = bnp.as_sequence_array(["acgttgta", "gcttca", "gttattc"], encoding=bnp.encodings.ACTGEncoding)

matrix = np.log([[0.1, 0.2],
                 [0.2, 0.3],
                 [0.4, 0.1],
                 [0.3, 0.4]])
pwm = bionumpy.sequence.position_weight_matrix.PositionWeightMatrix(matrix)
pwm("ac")
pwm(["ac", "cg"])
pwm.rolling_window(sequences)
pwm.rolling_window(sequences).max(axis=-1)


bnp.KmerEncoder(3).rolling_window(sequences)

fastq_entries_stream = bnp.open("example_data/reads.fq")
counts = np.zeros(4**3, dtype=int)
kmer_encoding = bnp.KmerEncoder(3)
for fastq_entries in fastq_entries_stream:
    kmer_hashes = kmer_encoding.rolling_window(fastq_entries.sequence)
    counts += np.bincount(kmer_hashes.ravel(), minlength=4**3)
counts
