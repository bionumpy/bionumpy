import numpy as np
from typing import Dict, List, Union
from .genomic_track import GenomicData
from .genome_context import GenomeContext
from ..datatypes import Interval
from ..io.indexed_fasta import IndexedFasta
from ..sequence import get_reverse_complement
from ..encodings import DNAEncoding, ACGTnEncoding
from ..encoded_array import as_encoded_array, EncodedArray, EncodedRaggedArray


def dna_encode(output):
    return as_encoded_array(output, ACGTnEncoding)


class GenomicSequence(GenomicData):
    '''Class to hold a genomic sequence as a GenomicArray'''
    def __init__(self, indexed_fasta: IndexedFasta, genome_context=None):
        self._genome_context = genome_context
        self._fasta = indexed_fasta

    @property
    def genome_context(self):
        if self._genome_context is None:
            return GenomeContext(self._fasta.get_contig_lengths())
        else:
            return self._genome_context

    def __repr__(self):
        return f'GenomicSequence over chromosomes: {list(self._fasta.keys())}'

    @classmethod
    def from_indexed_fasta(cls, indexed_fasta: IndexedFasta, genome_context=None) -> 'GenomicSequenceIndexedFasta':
        return GenomicSequenceIndexedFasta(indexed_fasta, genome_context)

    @classmethod
    def from_dict(cls, sequence_dict: Dict[str, str]):
        return GenomicSequenceDict(sequence_dict)

    def extract_chromsome(self, chromosome: Union[str, List[str]]) -> Union[EncodedArray, EncodedRaggedArray]:
        return dna_encode(self._fasta[chromosome])

    def _extract_intervals(self, intervals: Interval):
        return NotImplemented

    def _index_boolean(self, boolean_array: GenomicData):
        return self.extract_intervals(boolean_array.get_data(), stranded=False).ravel()

    def extract_intervals(self, intervals: Interval, stranded: bool = False) -> EncodedRaggedArray:
        sequences = self._extract_intervals(intervals)
        sequences = dna_encode(sequences)
        if stranded:
            sequences = np.where((intervals.strand == '+')[:, np.newaxis],
                                 sequences,
                                 get_reverse_complement(sequences))
        return sequences


class GenomicSequenceIndexedFasta(GenomicSequence):
    '''GenomicSeqeuence with an Indexed Fasta file as backend'''
    def _extract_intervals(self, intervals):
        return self._fasta.get_interval_sequences(intervals)


class GenomicSequenceDict(GenomicSequence):
    '''GenomicSeqeuence with a sequence dict as backend'''
    def __init__(self, sequence_dict: Dict[str, str]):
        self._dict = {name: as_encoded_array(sequence, target_encoding=ACGTnEncoding)
                      for name, sequence in sequence_dict.items()}

    def _extract_intervals(self, intervals):
        return as_encoded_array(
            [self._dict[interval.chromosome.to_string()][int(interval.start):int(interval.stop)]
             for interval in intervals],
            target_encoding=DNAEncoding)

    @property
    def genome_context(self):
        return GenomeContext({name: len(sequence) for name, sequence in self._dict.items()})