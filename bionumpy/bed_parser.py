from .npdataclass import NpDataClass, VarLenArray, SeqArray, npdataclass
from .parser import FileBuffer, NEWLINE
from .chromosome_provider import *
from dataclasses import dataclass
import numpy as np
        
class StrandEncoding:
    MIN_CODE = ord("+")
    @classmethod
    def from_bytes(cls, bytes_array):
        return (bytes_array & np.uint8(2)) >> np.uint8(1)

    @classmethod
    def to_bytes(cls, strands):
        return 2*strands + cls.MIN_CODE

class DigitEncoding:
    MIN_CODE = ord("0")
    @classmethod
    def from_bytes(cls, bytes_array):
        return bytes_array-cls.MIN_CODE

    @classmethod
    def to_bytes(cls, digits):
        return digits+cls.MIN_CODE

@npdataclass
class Interval:
    chromosome: SeqArray
    start: np.ndarray
    end: np.ndarray

    def in_interval(self, position):
        return (self.start <= position) & (position < self.end)

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

@npdataclass
class Variant:
    chromosome: SeqArray
    position: np.ndarray
    ref_seq: SeqArray
    alt_seq: SeqArray

class SNP(Variant):
    pass

class GenotypeEncoding:
    @classmethod
    def from_bytes(cls, bytes_array):
        assert bytes_array.shape[-1]==3
        return (bytes_array[..., 0]==ord("1"))+(bytes_array[..., 2]==ord("1")).astype("int")
