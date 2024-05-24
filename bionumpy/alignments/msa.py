from bionumpy import SequenceEntry
import numpy as np

class MultipleSequenceAlignment:
    def __init__(self, matrix, sequence_names):
        self.matrix = matrix
        self.sequence_names = sequence_names

    def matrix(self):
        return self.matrix

    @classmethod
    def from_sequence_entries(cls, entries: SequenceEntry):
        sequences = entries.sequence
        L = len(sequences[0])
        assert np.all(sequences.lengths == L)
        matrix = sequences.ravel().reshape(len(sequences), L)
        return cls(matrix, entries.name)

    def mask(self):
        return self.matrix != '-'