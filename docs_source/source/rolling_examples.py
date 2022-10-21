from bionumpy.rollable import RollableFunction
from bionumpy.sequences import as_encoded_array
import numpy as np


class MatchSequence(RollableFunction):
    def __init__(self, matching_sequence):
        self._matching_sequence = as_encoded_array(matching_sequence)
        self.window_size = self._matching_sequence.size
    def __call__(self, sequence):
        sequence = as_encoded_array(sequence)
        element_matches = sequence == self._matching_sequence
        return np.all(element_matches, axis=-1)


match = MatchSequence("CGGT")
match("CGGT")
match("CGGTA")
match.rolling_window("CGGTA")
match.rolling_window(["CGGTA", "ACGGTG"])
