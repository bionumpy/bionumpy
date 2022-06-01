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
    
    def __str__(self):
        strings = ("".join(chr(i) for i in array[:20]+"..."*(len(array)>20)) for array in self.tolist())
        seqs = ', '.join(seq for seq, _  in zip(strings, range(20)))
        trail = " ..." if len(self)>20 else ""
        return f"Sequences({seqs}{trail})"


def as_sequence_array(s, encoding=BaseEncoding):
    if isinstance(s, Sequence, Sequences):
        assert s.encoding==encoding
        return s
    elif isinstance(s, str):
        return Sequence.from_string(s)
    elif isinstance(s, list):
        return Sequences.from_strings(s)
