import os
import numpy as np
from typing import Dict
from pathlib import PurePath
from ..io import bnp_open, Bed6Buffer, BedBuffer, open_indexed
from ..io.files import buffer_types, BamBuffer, BamIntervalBuffer
from ..io.indexed_files import IndexBuffer, create_index
from .genomic_track import GenomicArray
from .genomic_intervals import GenomicIntervals, GenomicLocation, LocationEntry
from .genomic_sequence import GenomicSequence
from .annotation import GenomicAnnotation
from .genome_context import GenomeContext, keep_all, ignore_underscores
from ..encoded_array import as_encoded_array
from ..bnpdataclass import replace, BNPDataClass
from ..datatypes import BedGraph, Interval
from ..util.formating import table


class Genome:
    '''
    The entry point for working with genomic data. Supports reading data directly from
    file and into GenomicData objects, or creating GenomicData objects from `BNPDataClass` objects, 
    such `Interval`, `BedGraph` or `LocationEntry`.

    Supports both in memory or streamed/lazy data, controlled with the `stream` keyword. For streamed
    data it requires that the data in the file uses the same chromosome ordering as the genome. This
    will usually be the same as either what is provided in the `chrom.sizes` file, or a simple alphabetic
    ordering. If sort order errors occur, try using `sort_names=True` when creating the `Genome` object.
    '''
    def __init__(self, chrom_sizes: Dict[str, int], fasta_filename: str = None, sort_names=False, filter_function=lambda x: True):
        if isinstance(chrom_sizes, GenomeContext):
            self._genome_context = chrom_sizes
        else:
            if sort_names:
                chrom_sizes = {key: chrom_sizes[key] for key in sorted(chrom_sizes.keys())}
            self._genome_context = GenomeContext.from_dict(chrom_sizes, filter_function)
        self._fasta_filename = fasta_filename

    def with_ignored_added(self, ignored):
        return self.__class__(self._genome_context.with_ignored_added(ignored), self._fasta_filename)

    @classmethod
    def from_dict(cls, chrom_sizes: Dict[str, int], *args, **kwargs):
        return cls(chrom_sizes, *args, **kwargs)

    @classmethod
    def from_file(cls, filename: str, sort_names: bool=False, filter_function=ignore_underscores) -> 'Genome':
        """Read genome information from a 'chrom.sizes' or 'fa.fai' file

        File should contain rows of names and lengths of chromosomes.
        If a fasta file is provided, a fasta index will be written to 'fa.fai'
        if not already theree. The sequence will then be available through, 'genome.read_sequence()'
        without any 'filename' parameter.

        Parameters
        ----------
        cls :
        filename : str
        sort_names : bool Whether or not to sort the chromosome names

        Returns
        -------
        'Genome'
            A `Genome` object created from the read chromosome sizes

        """
        path = PurePath(filename)
        suffix = path.suffixes[-1]
        index_file_name = path.with_suffix(path.suffix + ".fai")
        fasta_filename = None
        if suffix in (".fa", ".fasta"):
            if not os.path.isfile(index_file_name):
                bnp_open(index_file_name, "w", buffer_type=IndexBuffer).write(create_index(path))
            fasta_filename = filename
            filename = index_file_name

        split_lines = (line.split()[:2] for line in open(filename))
        return cls({name: int(length) for name, length in split_lines}, fasta_filename=fasta_filename, sort_names=sort_names, filter_function=filter_function)

    @staticmethod
    def _open(filename, stream, buffer_type=None):
        f = bnp_open(filename, buffer_type=buffer_type)
        if stream:
            content = f.read_chunks()
        else:
            content = f.read()
        return content

    def get_track(self, bedgraph: BedGraph) -> GenomicArray:
        """Create a `GenomicArray` for the data in the bedgraph

        Parameters
        ----------
        bedgraph : BedGraph

        Returns
        -------
        GenomicArray
        """
        bedgraph = self._mask_data_on_extra_chromosomes(bedgraph)
        return GenomicArray.from_bedgraph(bedgraph, self._genome_context)

    def read_track(self, filename: str, stream: bool = False) -> GenomicArray:
        """Read a bedgraph from file and convert it to a `GenomicArray`

        If `stream` is `True` then read the bedgraph in chunks and get
        a lazy evaluated GenomicArray

        Parameters
        ----------
        filename : str
            Filename for the bedgraph file
        stream : bool
            Whether or not to read as stream

        """
        content = self._open(filename, stream)
        return self.get_track(content)

    def get_intervals(self, intervals: Interval, stranded: bool = False) -> GenomicIntervals:
        """Get genomic intervals from interval data. 

        Parameters
        ----------
        intervals : Interval
            Interval's or stream of Intervals's
        stranded : bool
            Wheter or not the intervals are stranded

        Returns
        -------
        GenomicIntervals
        """
        return GenomicIntervals.from_intervals(intervals, self._genome_context, is_stranded=stranded)

    def read_intervals(self, filename: str, stranded: bool = False, stream: bool = False, buffer_type=None) -> GenomicIntervals:
        """Read a bed file and represent it as `GenomicIntervals`

        If `stream` is `True` then read the bedgraph in chunks and get
        a lazy evaluated GenomicTrack

        Parameters
        ----------
        filename : str
            Filename for the bed file
        stranded : bool
            Whether or not to treat the intervals as stranded
        stream : bool
            Wheter to read as a stream

        Returns
        -------
        GenomicIntervals
        """
        path = PurePath(filename)
        suffix = path.suffixes[-1]
        if suffix == ".gz":
            suffix = path.suffixes[-2]
        if buffer_type is None:
            buffer_type = buffer_types[suffix]
            if buffer_type == BedBuffer and stranded: 
                buffer_type = Bed6Buffer
            if buffer_type == BamBuffer:
                buffer_type = BamIntervalBuffer
        content = self._open(filename, stream, buffer_type=buffer_type)
        return self.get_intervals(content, stranded)
    # return GenomicIntervals.from_intervals(content, self._chrom_sizes)

    def read_locations(self, filename: str, stranded: bool = False, stream: bool = False, has_numeric_chromosomes=False, buffer_type=None) -> GenomicLocation:
        """Read a set of locations from file and convert them to `GenomicLocation`        

        Locations can be treated as stranded or not stranded. If `has_numeric_chromosomes` then 'chr'
        will be prepended to all the chromosome names in order to make them compatible with the
        genome's chromosome names

        Parameters
        ----------
        filename : str
            Filename of the locations
        stranded : bool
            Whether or not the locations are stranded
        stream : bool
            Whether the data should be read as a stream
        has_numeric_chromosomes : bool
            Whether the file has the chromosomes as numbers

        Returns
        -------
        GenomicLocation
            (Stranded) GenomicLocation

        """
        assert not (stream and has_numeric_chromosomes)
        f = bnp_open(filename, buffer_type=buffer_type)
        data = f.read_chunks()
        if not stream:
            data_list = list(data)
            if len(data_list) == 0:
                data = LocationEntry.empty()
            else:
                data = np.concatenate(data_list)

        return self.get_locations(data, has_numeric_chromosomes=has_numeric_chromosomes)

    def _mask_data_on_extra_chromosomes(self, data, chromosome_field_name='chromosome'):
        if (not isinstance(data, BNPDataClass)) or len(data)==0:
            return data
        encoded_chromosomes = self._genome_context.encoding.encode(getattr(data, chromosome_field_name))
        data = replace(data, **{chromosome_field_name: encoded_chromosomes})
        mask = self._genome_context.is_included(encoded_chromosomes)
        return data[mask]

    def get_locations(self, data: LocationEntry, has_numeric_chromosomes=False) -> GenomicLocation:
        """Create `GenomicLocation`s from location entries

        Parameters
        ----------
        data : LocationEntry
            BNPDataClass of location entries

        Returns
        -------
        GenomicLocation

        """
        if has_numeric_chromosomes:
            data = replace(data, chromosome=as_encoded_array(['chr'+chromosome.to_string() for chromosome in data.chromosome]))
        data = self._mask_data_on_extra_chromosomes(data)
        return GenomicLocation.from_data(data, self._genome_context)

    def read_sequence(self, filename: str = None) -> GenomicSequence:
        """Read the genomic sequence from file.

        If a `fasta` file was used to create the Genome object, `filename` can be `None`
        in which case the sequence will be read from that fasta file

        Parameters
        ----------
        filename : str

        Returns
        -------
        GenomicSequence

        """
        
        if filename is None:
            assert self._fasta_filename is not None
            filename = self._fasta_filename
        return GenomicSequence.from_indexed_fasta(open_indexed(filename), genome_context=self._genome_context)

    def read_annotation(self, filename: str) -> GenomicAnnotation:
        """Read genomic annotation from a 'gtf' or 'gff' file

        The annotation file should have 'gene', 'transcript' and 'exon' entries.

        Parameters
        ----------
        filename : str
            Filename of the annotation file


        Returns
        -------
        GenomicAnnotation
            `GenomicAnnotation` containing the genes, transcripts and exons

        """
        
        gtf_entries = self._open(filename, stream=False)
        return GenomicAnnotation.from_gtf_entries(gtf_entries, self._genome_context)

    def __repr__(self):
        return f"{self.__class__.__name__}(" + repr(self._genome_context) + ")"

    def __str__(self):
        return table(
            ((key, value) for key, value in self._genome_context.chrom_sizes.items() if '_' not in key),
            headers=["Chromosome", "Size"])

    def get_genome_context(self):
        return self._genome_context

    @property
    def size(self):
        '''The size of the genome'''
        return self._genome_context.size
