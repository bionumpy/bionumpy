from bionumpy import as_sequence_array
import numpy as np

sequences = as_sequence_array([
    "acttctagtcggatctgt",
    "cgatcggatcgttgc",
    "ctattcgattcagtgtggatgc",
    "cttctagtat"])
trimmed_sequences = sequences[:, 3:-3]
g_mask = sequences == as_sequence_array("g")
np.sum(g_mask, axis=-1)
# [4 5 6 1]
