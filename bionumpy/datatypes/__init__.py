from ..typing import SequenceID
import numpy as np
from typing import List, Union, Optional
from ..encodings import (CigarOpEncoding, BamEncoding, QualityEncoding,
                         CigarEncoding, StrandEncoding)
from ..encodings.vcf_encoding import PhasedGenotypeRowEncoding, GenotypeRowEncoding, PhasedHaplotypeRowEncoding
from ..bnpdataclass import bnpdataclass, BNPDataClass
from .gtf import GFFEntry, GFFExonEntry, GFFGeneEntry, GFFTranscriptEntry, GTFEntry
from ..config import STRING_ARRAY

if not STRING_ARRAY:
    SequenceID = str


@bnpdataclass
class LocationEntry:
    chromosome: SequenceID
    position: int


@bnpdataclass
class StrandedLocationEntry(LocationEntry):
    strand: StrandEncoding


@bnpdataclass
class BedGraph:
    chromosome: SequenceID
    start: int
    stop: int
    value: float


@bnpdataclass
class RawSeqeuence:
    sequence: str


@bnpdataclass
class SequenceEntry:
    name: SequenceID
    sequence: str


@bnpdataclass
class SequenceEntryWithQuality(SequenceEntry):
    quality: QualityEncoding


@bnpdataclass
class Interval:
    chromosome: SequenceID
    start: int
    stop: int


@bnpdataclass
class StrandedInterval(Interval):
    strand: StrandEncoding


@bnpdataclass
class NamedInterval(Interval):
    name: SequenceID


@bnpdataclass
class Bed6(NamedInterval):
    score: Optional[int]
    strand: StrandEncoding


@bnpdataclass
class NarrowPeak(Bed6):
    signal_value: float
    p_value: float
    q_value: float
    summit: int


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
    chromosome: SequenceID
    position: int
    ref_seq: str
    alt_seq: str

    def is_snp(self):
        return (self.ref_seq.lengths == 1) & (self.alt_seq.lengths == 1)


@bnpdataclass
class VCFEntry:
    chromosome: SequenceID
    position: int
    id: str
    ref_seq: str
    alt_seq: str
    quality: str
    filter: str
    info: Union[BNPDataClass, str]
    # genotypes: str

    def is_snp(self):
        return (self.ref_seq.lengths == 1) & (self.alt_seq.lengths == 1)


@bnpdataclass
class VCFWithInfoAsStringEntry:
    chromosome: SequenceID
    position: int
    id: str
    ref_seq: str
    alt_seq: str
    quality: str
    filter: str
    info: str
    # genotypes: str

    def is_snp(self):
        return (self.ref_seq.lengths == 1) & (self.alt_seq.lengths == 1)


@bnpdataclass
class VCFEntryWithGenotypes(VCFEntry):
    genotype: List[str]


@bnpdataclass
class VCFGenotypeEntry(VCFEntry):
    genotypes: GenotypeRowEncoding


@bnpdataclass
class PhasedVCFGenotypeEntry(VCFEntry):
    genotypes: PhasedGenotypeRowEncoding


@bnpdataclass
class PhasedVCFHaplotypeEntry(VCFEntry):
    genotypes: PhasedHaplotypeRowEncoding


class SNP(Variant):
    pass


@bnpdataclass
class SAMEntry:
    name: SequenceID
    flag: int
    chromosome: SequenceID
    position: int
    mapq: int
    cigar: str
    next_chromosome: str
    next_position: int
    length: int
    sequence: str
    quality: str
    extra: str

@bnpdataclass
class BamEntry:
    chromosome: SequenceID
    name: SequenceID
    flag: int
    position: int
    mapq: int
    cigar_op: CigarOpEncoding
    cigar_length: CigarEncoding
    sequence: BamEncoding
    quality: QualityEncoding


@bnpdataclass
class ChromosomeSize:
    name: str
    size: int


@bnpdataclass
class GfaPath:
    name: str
    node_ids: List[int]
    directions: List[int]


@bnpdataclass
class PairsEntry:
    """https://pairtools.readthedocs.io/en/latest/formats.html"""
    read_id: str
    chrom1: SequenceID
    pos1: int
    chrom2: SequenceID
    pos2: int
    strand1: StrandEncoding
    strand2: StrandEncoding

