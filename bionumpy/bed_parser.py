from .npdataclass import NpDataClass, VarLenArray, SeqArray
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

class Interval:
    def __init__(self, data):
        self.data = np.asanyarray(data)
        assert self.data.shape[-1] == 2
        
        self.start = self.data[..., 0]
        self.end = self.data[..., 1]

    def __repr__(self):
        return repr(self.data)

    def __eq__(self, other):
        return np.all(self.data==other.data)

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

"""
class FullBedFile(ChromosomeDictProvider):
    def __init__(self, chrom_dict, all_intervals):
        self._chrom_dict = chrom_dict
        self._all_intervals = all_intervals
        
    def __getitem__(self, chrom):
        start, end = chrom_dict[chrom]
        return all_intervals[start:end]

    @classmethod
    def from_bed_buffer_stream(cls, bed_buffer_stream):
        cur_offset = 0
        all_intervals = []
        chrom_starts = {}
        chrom_ends = {}
        for bed_buffer in bed_buffer_stream:
            chromosomes = cur_buffer.get_chromosomes()
            if not len(chromosomes):
                break
            if not cls._is_same_chromosome(last_chromosome, chromosomes[0]):
                chrom_starts[cls.get_chrom_name(chromosomes[0])] == cur_offset
            chromosome_changes = cls._get_chromosome_changes(chromosomes)
            chrom_starts.update((cls.get_chrom_name(chromosomes[i]), i+cur_offset)
                                for i in chromosome_changes)
            all_intervals.append(cur_buffer.get_intervals())
            last_chromosome = chromosomes[-1]
            cur_offset += len(chromosomes)
        chrom_ends[cls.get_chrom_name(last_chromosome)] == cur_offset
        interval_dict = {chrom: (chrom_starts[chrom], chrom_ends[chrom]) for chrom in chrom_starts.keys()}
        return cls(interval_dict, np.concatenate(intervals))
"""

@dataclass
class Variant(NpDataClass):
    chromosome: SeqArray
    position: np.ndarray
    ref_seq: SeqArray
    alt_seq: SeqArray

    def __eq__(self, other):
        return all(np.all(np.equal(s, o)) for s, o in zip(self.shallow_tuple(), other.shallow_tuple()))

class SNP(Variant):
    pass

class GenotypeEncoding:
    @classmethod
    def from_bytes(cls, bytes_array):
        assert bytes_array.shape[-1]==3
        return (bytes_array[..., 0]==ord("1"))+(bytes_array[..., 2]==ord("1")).astype("int")
