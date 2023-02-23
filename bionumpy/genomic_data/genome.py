import os
import numpy as np
from typing import Dict
from pathlib import PurePath
from ..io import bnp_open, Bed6Buffer, BedBuffer, open_indexed
from ..io.files import buffer_types, BamBuffer, BamIntervalBuffer
from ..io.indexed_files import IndexBuffer, create_index
from .genomic_track import GenomicArray
from .genomic_intervals import GenomicIntervals, GenomicLocation
from .genomic_sequence import GenomicSequence
from .annotation import GenomicAnnotation
from .geometry import Geometry, StreamedGeometry
from ..encodings.string_encodings import StringEncoding
from ..util.formating import table
from ..encoded_array import as_encoded_array
from ..bnpdataclass import replace, BNPDataClass


class Genome:
    '''Should return GenomicIntervals, GenomicTrack, GenomicMask'''
    def __init__(self, chrom_sizes: Dict[str, int], fasta_filename: str = None, sort_names=False):
        if sort_names:
            chrom_sizes = {key: chrom_sizes[key] for key in sorted(chrom_sizes.keys())}
        self._chrom_sizes = chrom_sizes
        self._geometry = Geometry(self._chrom_sizes)
        self._chromosome_encoding = StringEncoding(as_encoded_array(list(chrom_sizes.keys())))
        self._streamed_geometry = StreamedGeometry(self._chrom_sizes)
        self._fasta_filename = fasta_filename
        self._include_mask = np.array(['_' not in name for name in chrom_sizes])

    @classmethod
    def from_file(cls, filename: str, sort_names=False) -> 'Genome':
        """Read genome information from a 'chrom.sizes' or 'fa.fai' file

        File should contain rows of names and lengths of chromosomes

        Parameters
        ----------
        cls :
        filename : str

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
        return cls({name: int(length) for name, length in split_lines}, fasta_filename=fasta_filename, sort_names=sort_names)

    @staticmethod
    def _open(filename, stream, buffer_type=None):
        f = bnp_open(filename, buffer_type=buffer_type)
        if stream:
            content = f.read_chunks()
        else:
            content = f.read()
        return content

    def get_track(self, bedgraph):
        bedgraph = self._mask_data_on_extra_chromosomes(bedgraph)
        return GenomicArray.from_bedgraph(bedgraph, self._chrom_sizes)

    def read_track(self, filename: str, stream: bool=False) -> GenomicArray:
        """Read a bedgraph from file and convert it to a `GenomicTrack`

        If `stream` is `True` then read the bedgraph in chunks and get
        a lazy evaluated GenomicTrack

        Parameters
        ----------
        filename : str
            Filename for the bedgraph file
        stream : bool
            Whether or not to read as stream

        """
        content = self._open(filename, stream)
        return self.get_track(content)
        # return GenomicArray.from_bedgraph(content, self._chrom_sizes)

    def __read_mask(self, filename: str , stream: bool = False) -> GenomicArray:
        """Read a bed file and convert it to a `GenomicMask` of areas covered by an interval

        If `stream` is `True` then read the bedgraph in chunks and get
        a lazy evaluated GenomicTrack

        Parameters
        ----------
        filename : str
            Filename for the bed file
        stream : bool
            Wheter to read as a stream

        """
        content = self._open(filename, stream)
        geom = self._streamed_geometry if stream else self._geometry
        return geom.get_mask(content)

    def get_intervals(self, intervals, stranded=False):
        intervals = self._mask_data_on_extra_chromosomes(intervals)
        return GenomicIntervals.from_intervals(intervals, self._chrom_sizes, is_stranded=stranded)

    def read_intervals(self, filename: str, stranded: bool = False, stream: bool = False) -> GenomicIntervals:
        """Read a bed file and represent it as `GenomicIntervals`

        If `stream` is `True` then read the bedgraph in chunks and get
        a lazy evaluated GenomicTrack

        Parameters
        ----------
        filename : str
            Filename for the bed file
        stream : bool
            Wheter to read as a stream

        """
        path = PurePath(filename)
        suffix = path.suffixes[-1]
        if suffix == ".gz":
            suffix = path.suffixes[-2]
        buffer_type = buffer_types[suffix]
        if buffer_type == BedBuffer and stranded: 
            buffer_type = Bed6Buffer
        if buffer_type == BamBuffer:
            buffer_type = BamIntervalBuffer
        content = self._open(filename, stream, buffer_type=buffer_type)
        return self.get_intervals(content, stranded)
    # return GenomicIntervals.from_intervals(content, self._chrom_sizes)

    def _mask_data_on_extra_chromosomes(self, data, chromosome_field_name='chromosome'):
        if not isinstance(data, BNPDataClass):
            return data
        encoded_chromosomes = self._chromosome_encoding.encode(getattr(data, chromosome_field_name))
        data = replace(data, **{chromosome_field_name: encoded_chromosomes})
        mask = self._include_mask[encoded_chromosomes.raw()]
        return data[mask]

    def get_locations(self, data):
        data = self._mask_data_on_extra_chromosomes(data)
        return GenomicLocation.from_data(data, self._chrom_sizes)

    def read_sequence(self, filename: str = None) -> GenomicSequence:
        if filename is None:
            assert self._fasta_filename is not None
            filename = self._fasta_filename 
        return GenomicSequence.from_indexed_fasta(open_indexed(filename))

    def read_annotation(self, filename: str) -> GenomicAnnotation:
        gtf_entries = self._open(filename, stream=False)
        return GenomicAnnotation.from_gtf_entries(gtf_entries, self._chrom_sizes)

    def __repr__(self):
        return f"{self.__class__.__name__}(" + repr(self._chrom_sizes) + ")"

    def __str__(self):
        return table(
            ((key, value) for key, value in self._chrom_sizes.items() if '_' not in key),
            headers=["Chromosome", "Size"])

    @property
    def size(self):
        return sum(self._chrom_sizes.values())
