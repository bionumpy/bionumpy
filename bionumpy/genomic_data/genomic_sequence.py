import numpy as np
from typing import Dict
from .genomic_track import GenomicData, GenomeContext
from ..io.indexed_fasta import IndexedFasta
from ..sequence import get_reverse_complement
from ..encodings import DNAEncoding, ACGTnEncoding
from ..encoded_array import as_encoded_array


def dna_encode(output):
    return as_encoded_array(output, ACGTnEncoding)


class GenomicSequence(GenomicData):
    def __init__(self, indexed_fasta: IndexedFasta):
        self._fasta = indexed_fasta

    @property
    def genome_context(self):
        return GenomeContext(self._fasta.get_contig_lengths())

    def __repr__(self):
        return f'GenomicSequence over chromosomes: {list(self._fasta.keys())}'

    @classmethod
    def from_indexed_fasta(cls, indexed_fasta: IndexedFasta):
        return GenomicSequenceIndexedFasta(indexed_fasta)

    @classmethod
    def from_dict(cls, sequence_dict: Dict[str, str]):
        return GenomicSequenceDict(sequence_dict)

    def extract_chromsome(self, chromosome):
        return dna_encode(self._fasta[chromosome])

    def _extract_intervals(self, intervals):
        return NotImplemented

    def _index_boolean(self, boolean_array):
        return self.extract_intervals(boolean_array.get_data(), stranded=False).ravel()

    def extract_intervals(self, intervals, stranded: bool = False):
        sequences = self._extract_intervals(intervals)
        if stranded:
            sequences = np.where((intervals.strand == '+')[:, np.newaxis],
                                 sequences,
                                 get_reverse_complement(sequences))
        return dna_encode(sequences)


class GenomicSequenceIndexedFasta(GenomicSequence):

    def _extract_intervals(self, intervals):
        return self._fasta.get_interval_sequences(intervals)


class GenomicSequenceDict(GenomicSequence):
    def __init__(self, sequence_dict: Dict[str, str]):
        self._dict = {name: as_encoded_array(sequence, target_encoding=ACGTnEncoding)
                      for name, sequence in sequence_dict.items()}

    def _extract_intervals(self, intervals):
        return as_encoded_array(
            [self._dict[interval.chromosome.to_string()][int(interval.start):int(interval.stop)]
             for interval in intervals],
            target_encoding=DNAEncoding)