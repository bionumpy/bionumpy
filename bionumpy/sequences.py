import numpy as np
from .encodings import ACTGTwoBitEncoding, BaseEncoding
from npstructures import RaggedArray


class Sequences(RaggedArray):
    def __init__(self, data, shape=None, encoding=BaseEncoding):
        super().__init__(data, shape, dtype=np.uint8)
        self.encoding = encoding

    @classmethod
    def from_sequences(cls, sequences):
        return cls([[ord(c) for c in seq] for seq in sequences])

    def to_sequences(self):
        return ["".join(chr(i) for i in array) for array in self.tolist()]


class _Sequences:
    def __init__(self, sequences, intervals, encoding=BaseEncoding):
        assert intervals[0].size == intervals[1].size
        self.sequences = sequences
        self.intervals = intervals
        self.intervals_start = intervals[0]
        self.intervals_end = intervals[1]
        self.encoding = encoding

    def __len__(self):
        return len(self.intervals_start)

    def __getitem__(self, i):
        start, end = (self.intervals_start[i], self.intervals_end[i])
        return self.sequences[start:end]

    def __iter__(self):
        return (self.sequences[start:end] for start, end in zip(*self.intervals))

    def __repr__(self):
        return "Seqs(%s, %s)" % (str(self.sequences), self.intervals)
