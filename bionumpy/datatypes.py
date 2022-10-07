import numpy as np
from typing import List
from .encodings import QualityEncoding
from .encodings.base_encoding import CigarEncoding
from .encodings.alphabet_encoding import CigarOpArray, BamArray
from .bnpdataclass import bnpdataclass


@bnpdataclass
class RawSeqeuence:
    sequence: str


@bnpdataclass
class SequenceEntry:
    name: str
    sequence: str


@bnpdataclass
class SequenceEntryWithQuality:
    name: str
    sequence: str
    quality: QualityEncoding


@bnpdataclass
class Interval:
    chromosome: str
    start: int
    end: int


@bnpdataclass
class Bed6(Interval):
    name: str
    score: int
    strand: str


@bnpdataclass
class Bed12(Bed6):
    thick_start: int
    thick_end: int
    item_rgb: str
    block_count: int
    block_sizes: List[int]
    block_starts: List[int]


@bnpdataclass
class Variant:
    chromosome: str
    position: int
    ref_seq: str
    alt_seq: str


@bnpdataclass
class VCFEntry:
    chromosome: str
    position: int(-1)
    id: str
    ref_seq: str
    alt_seq: str
    # quality: int
    # filter: str
    # info: str


@bnpdataclass
class VariantWithGenotypes(Variant):
    genotypes: int


class SNP(Variant):
    pass


class SortedIntervals:
    def __init__(self, data):
        self.data = np.asanyarray(data)
        assert data.shape[-1] == 2
        assert len(data.shape) == 2

        self.starts = self.data[..., 0]
        self.ends = self.data[..., 1]

    def in_intervals(self, position):
        idx = np.minimum(
            np.searchsorted(self.starts, position, side="left"), self.starts.size - 1
        )
        return (position >= self.starts[idx]) & (position < self.ends[idx])

    @classmethod
    def concatenate(cls, elements):
        return cls(np.vstack([element.data for element in elements]))


@bnpdataclass
class GFFEntry:
    chromosome: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: int
    atributes: str


@bnpdataclass
class SAMEntry:
    name: str
    flag: int
    chromosome: str
    position: int
    mapq: int
    cigar: str
    next_chromosome: str
    next_position: int
    length: int
    sequence: str
    quality: str


@bnpdataclass
class BamEntry:
    chromosome: str
    name: str
    flag: int
    position: int
    mapq: int
    cigar_op: CigarOpArray
    cigar_length: CigarEncoding
    sequence: BamArray
    quality: QualityEncoding


@bnpdataclass
class ChromosomeSize:
    name: str
    size: int
