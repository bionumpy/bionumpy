from .encodings import ACTGTwoBitEncoding, BaseEncoding

class Sequences:
    def __init__(self, sequences, intervals, encoding=BaseEncoding):
        assert intervals[0].size == intervals[1].size
        self.sequences = sequences
        self.intervals = intervals
        self.intervals_start = intervals[0]
        self.intervals_end = intervals[1]
        self.encoding = encoding

    def __len__(self):
        return len(self.offsets)+1

    def __getitem__(self, i):
        start, end = (0, len(self.sequences))
        if i > 0:
            start = self.offsets[i-1]
        if i < len(self.offsets):
            end = self.offsets[i]
        return self.sequences[start:end]

    def __repr__(self):
        return "Seqs(%s, %s)" % (str(self.sequences), self.intervals)

