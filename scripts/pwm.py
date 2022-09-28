import bionumpy as bnp
import matplotlib.pyplot as plt
from bionumpy.jaspar import read_jaspar_matrix
from bionumpy.position_weight_matrix import PositionWeightMatrix, pwm_from_counts
from bionumpy.encodings.alphabet_encoding import get_alphabet_array_class

alphabet, matrix = read_jaspar_matrix("example_data/MA0080.1.jaspar")
pwm = pwm_from_counts(matrix)
arrayclass = get_alphabet_array_class(alphabet)
motif_score = PositionWeightMatrix(pwm, arrayclass)

entries = bnp.open("example_data/big.fq.gz").read()
print(entries)
print(repr(entries.sequence))
scores = motif_score.rolling_window(entries.sequence)
plt.hist(scores.max(axis=-1));plt.show()
#print(scores.max(axis=-1))
