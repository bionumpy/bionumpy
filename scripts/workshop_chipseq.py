import bionumpy as bnp
from pyjaspar import jaspardb
from bionumpy.sequence.position_weight_matrix import PWM, get_motif_scores
import numpy as np
jdb_obj = jaspardb(release="JASPAR2020")
human_motifs = jdb_obj.fetch_motifs(collection="CORE", species=["9606"])
peaks = bnp.open("/home/knut/Downloads/ENCFF843VHC.bed.gz").read_chunks()
reference_genome = bnp.open_indexed("/home/knut/Data/hg38.fa")
multistream = bnp.MultiStream(reference_genome.get_contig_lengths(),
                              intervals=peaks,
                              sequence=reference_genome)

sequence_stream = bnp.sequence.get_sequences(multistream.sequence,
                                             multistream.intervals)

sequences = np.concatenate(list(sequence_stream))

for motif in human_motifs:
    pwm = PWM.from_dict(motif.pwm)
    motif_scores = get_motif_scores(sequences, pwm)
    print(motif_scores)
