import bionumpy as bnp
from pyjaspar import jaspardb
from bionumpy.sequence.position_weight_matrix import PWM, get_motif_scores
import numpy as np
import plotly.express as plx
jdb_obj = jaspardb(release="JASPAR2020")
human_motifs = jdb_obj.fetch_motifs(collection="CORE", species=["9606"])
peaks = bnp.open("/home/knut/Downloads/ENCFF843VHC.bed.gz").read()
sorted_peaks = bnp.arithmetics.sort_all_intervals(peaks)

reference_genome = bnp.open_indexed("/home/knut/Data/hg38.fa")
contig_lenghts = bnp.MultiStream.sort_dict_by_key(reference_genome.get_contig_lengths())
multistream = bnp.MultiStream(contig_lenghts,
                              intervals=sorted_peaks,
                              sequence=reference_genome)

sequence_stream = bnp.sequence.get_sequences(multistream.sequence,
                                             multistream.intervals)

sequences = np.concatenate(list(sequence_stream))

counts = []
lengths = []
for i, motif in enumerate(human_motifs):
    if i % 10 == 0:
        print(i)
    pwm = PWM.from_dict(motif.pwm)
    motif_scores = get_motif_scores(sequences, pwm)
    adjusted_scores = motif_scores-np.logaddexp(motif_scores, motif.length*np.log(0.25))
    counts.append((adjusted_scores.max(axis=-1) > np.log(0.9)).sum())
    lengths.append(motif.length)

plx.scatter(lengths, counts)
names = [m.name for m in human_motifs]
args = np.argsort(counts)
sorted_names = [names[i] for i in args]
