from bionumpy import as_encoded_array
import numpy as np

sequences = as_encoded_array([
    "acttctagtcggatctgt",
    "cgatcggatcgttgc",
    "ctattcgattcagtgtggatgc",
    "cttctagtat"])
trimmed_sequences = sequences[:, 3:-3]
g_mask = sequences == "g"
np.sum(g_mask, axis=-1)
# [4 5 6 1]
