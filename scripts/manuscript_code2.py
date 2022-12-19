from bionumpy.sequence.position_weight_matrix import PWM, get_motif_scores
import bionumpy as bnp
import numpy as np
from pyjaspar import jaspardb
import plotly.express as plx

# Read peaks and reference genom
peaks = bnp.open("ctcf.bed.gz").read()
reference_genome = bnp.open_indexed("hg38.fa")

# Change start and end position of peaks (100 bp around centre)
midpoints = (peaks.start+peaks.stop)//2
peaks.start = midpoints - 50
peaks.stop = midpoints + 50

# Fetch sequences within each peak
peak_sequences = reference_genome.get_interval_sequences(peaks)

# Fetch a motif from Jaspar and read it into BioNumPy
jaspar_object = jaspardb(release="JASPAR2020")
ctcf_motif = jaspar_object.fetch_motifs_by_name('CTCF')[0]
ctcf_pwm = PWM.from_dict(ctcf_motif.pwm)

# Get motif scores and make a boolean mask of likely binding sites
motif_scores = get_motif_scores(peak_sequences, ctcf_pwm)
has_motif_match = motif_scores > np.log(4)

# Plot mean matches per base (axis=0)
fig = plx.line(y=np.mean(has_motif_match, axis=0), template='seaborn',
               labels={"x": "Position in peak",
                       "y": "Ratio of peaks with motif match"},
               )
fig.write_image("motif_matches.png")
