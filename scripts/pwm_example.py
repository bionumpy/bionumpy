import bionumpy as bnp
from bionumpy.jaspar import read_jaspar_matrix
from bionumpy.position_weight_matrix import PositionWeightMatrix, pwm_from_counts
from bionumpy.encodings.alphabet_encoding import AlphabetEncoding

# Read the alphabet and counts from jaspar file
alphabet, matrix = read_jaspar_matrix("example_data/MA0080.1.jaspar")

# Convert counts to position weight matrix
pwm = pwm_from_counts(matrix)

# Make an array-class for the alphabet
encoding = AlphabetEncoding(alphabet)

# Get the motif score function
motif_score = PositionWeightMatrix(pwm, encoding)

#Get reads
entries = bnp.open("example_data/big.fq.gz").read()

# Calculate the motif score for each valid window
scores = motif_score.rolling_window(entries.sequence)

# Get a histogram of the max-score for each read
bnp.histogram(scores.max(axis=-1))
