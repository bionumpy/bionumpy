import numpy as np

from .npdataclass import npdataclass, SeqArray

@npdataclass
class SequenceEntry:
    name: SeqArray
    sequence: SeqArray

@npdataclass
class SequenceEntryWithQuality:
    name: SeqArray
    sequence: SeqArray
    quality: SeqArray

@npdataclass
class Interval:
    chromosome: SeqArray
    start: np.ndarray
    end: np.ndarray

    def in_interval(self, position):
        return (self.start <= position) & (position < self.end)

@npdataclass
class Variant:
    chromosome: SeqArray
    position: np.ndarray
    ref_seq: SeqArray
    alt_seq: SeqArray

    def is_snp(self):
        return np.logical_and(self.ref_seq.shape.lengths == 1, self.alt_seq.shape.lengths == 1)

class SNP(Variant):
    pass


"""
class SortedIntervals:
    def __init__(self, data):
        self.data = np.asanyarray(data)
        assert data.shape[-1] == 2
        assert len(data.shape) == 2
        
        self.starts = self.data[..., 0]
        self.ends = self.data[..., 1]

    def in_intervals(self, position):
        idx = np.minimum(np.searchsorted(self.starts, position, side="left"), self.starts.size-1)
        return (position >=self.starts[idx]) & (position < self.ends[idx])

    @classmethod
    def concatenate(cls, elements):
        return cls(np.vstack([element.data for element in elements]))
"""
