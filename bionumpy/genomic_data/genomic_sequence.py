from .genomic_track import GenomicData
from ..io.indexed_fasta import IndexedFasta


class GenomicSequence(GenomicData):
    def __init__(self, indexed_fasta: IndexedFasta):
        self._fasta = indexed_fasta

    def __repr__(self):
        return f'GenomicSequence over chromosomes: {list(self._fasta.keys())}'

    @classmethod
    def from_indexed_fasta(cls, indexed_fasta: IndexedFasta):
        return cls(indexed_fasta)

    def extract_chromsome(self, chromosome):
        return self._fasta[chromosome]

    def extract_intervals(self, intervals):
        return self._fasta.get_interval_sequences(intervals)
