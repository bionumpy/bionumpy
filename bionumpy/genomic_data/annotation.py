from ..datatypes import GTFEntry
from .genomic_intervals import GenomicIntervalsFull
from .genome_context_base import GenomeContextBase


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
    ''' Class to hold a genomic annotations. Basically just a holder for gene, transcript and exon intervals.'''

    def __init__(self, data, genome_context):
        # self._genes, transcripts, exons, data):
        # self._genes = genes
        # self._transcripts = transcripts
        # self._exons = exons
        self._exons = None
        self._transcripts = None
        self._genes = None
        self._data = data
        self._genome_context = genome_context

    def __str__(self):
        return f'GenomicAnnotation with {len(self._genes)} genes, {len(self._transcripts)} transcripts and {len(self._exons)} exons'

    def __repr__(self):
        return f'GenomicAnnotation(genome_context={self._genome_context}, data={self._data})'

    @property
    def genes(self) -> Genes:
        if self._genes is None:
            self._genes = Genes(self._data.get_genes(), self._genome_context, True)
        return self._genes

    @property
    def transcripts(self) -> Transcripts:
        if self._transcripts is None:
            self._transcripts = Transcripts(self._data.get_transcripts(), self._genome_context, True)
        return self._transcripts

    @property
    def exons(self) -> Exons:
        if self._exons is None:
            self._exons = Exons(self._data.get_exons(), self._genome_context, True)
        return self._exons

    @classmethod
    def from_gtf_entries(cls, gtf_entries: GTFEntry, genome_context: GenomeContextBase) -> 'GenomicAnnotation':
        """Get a GenomicAnnotation from a set of gtf/gff entries

        Parameters
        ----------
        gtf_entries : GTFEntry
            BNPDataClass with gtf entries
        genome_context : GenomeContextBase
            The genome context

        Returns
        -------
        'GenomicAnnotation'
            Annotatin object holding the genes, transcripts and exons
        """
        return cls(gtf_entries, genome_context)# Genes(gtf_entries.get_genes(), genome_context, True),
    #Transcripts(gtf_entries.get_transcripts(), genome_context, True),
    #               Exons(gtf_entries.get_exons(), genome_context, True),
    #               gtf_entries)
