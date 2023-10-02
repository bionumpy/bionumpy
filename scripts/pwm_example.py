import bionumpy as bnp
from bionumpy.io.motifs import read_motif
from bionumpy.sequence.position_weight_matrix import PWM, get_motif_scores

# Read the alphabet and counts from jaspar file
pwm = read_motif("example_data/MA0080.1.jaspar")

# Convert counts to position weight matrix
# pwm = PWM.from_counts(motif)

# Get reads
entries = bnp.open("example_data/big.fq.gz").read()

# Calculate the motif score for each valid window
scores = get_motif_scores(entries.sequence, pwm)

# Get a histogram of the max-score for each read
bnp.histogram(scores.max(axis=-1))
