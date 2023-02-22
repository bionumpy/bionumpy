from typing import Dict
from ..datatypes import GTFEntry
from .genomic_intervals import GenomicIntervalsFull


class Genes(GenomicIntervalsFull):

    @property
    def gene_id(self):
        return self._intervals.gene_id


class Transcripts(Genes):
    @property
    def transcript_id(self):

        return self._intervals.transcript_id


class Exons(Transcripts):
    @property
    def exon_id(self):
        return self._intervals.exon_id


class GenomicAnnotation:

    @classmethod
    def __init__(self, genes, transcripts, exons, data):
        self._genes = genes
        self._transcripts = transcripts
        self._exons = exons
        self._data = data

    @property
    def genes(self):
        return self._genes

    @property
    def transcripts(self):
        return self._transcripts

    @property
    def exons(self):
        return self._exons

    @classmethod
    def from_gtf_entries(cls, gtf_entries: GTFEntry, chromosome_sizes: Dict[str, int]):
        return cls(Genes(gtf_entries.get_genes(), chromosome_sizes, True),
                   Transcripts(gtf_entries.get_transcripts(), chromosome_sizes, True),
                   Exons(gtf_entries.get_exons(), chromosome_sizes, True),
                   gtf_entries)
