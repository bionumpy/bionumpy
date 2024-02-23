import numpy as np
import bionumpy as bnp
from bionumpy.io.motifs import read_motif
from bionumpy.sequence.position_weight_matrix import get_motif_scores


def read_motif_scores(reads_filename: str, motif_filename: str) -> np.ndarray:
    # Read the alphabet and counts from jaspar file
    pwm = read_motif(motif_filename)

    # Get reads
    entries = bnp.open(reads_filename).read()

    # Calculate the motif score for each valid window
    scores = get_motif_scores(entries.sequence, pwm)

    # Get a histogram of the max-score for each read
    return bnp.histogram(scores.max(axis=-1))


if __name__ == "__main__":
    read_motif_scores("example_data/big.fq.gz", "example_data/MA0080.1.jaspar")