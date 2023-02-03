from typing import Dict
from ..io import bnp_open, Bed6Buffer, BedBuffer
from .genomic_track import GenomicTrack
from .genomic_intervals import GenomicIntervals
from .geometry import Geometry, StreamedGeometry


class Genome:
    '''Should return GenomicIntervals, GenomicTrack, GenomicMask'''
    def __init__(self, chrom_sizes: Dict[str, int]):
        self._chrom_sizes = chrom_sizes
        self._geometry = Geometry(self._chrom_sizes)
        self._streamed_geometry = StreamedGeometry(self._chrom_sizes)

    @classmethod
    def from_file(cls, filename: str) -> 'Genome':
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
        split_lines = (line.split()[:2] for line in open(filename))
        return cls({name: int(length) for name, length in split_lines})

    @staticmethod
    def _open(filename, stream, buffer_type=None):
        f = bnp_open(filename, buffer_type=buffer_type)
        if stream:
            content = f.read_chunks()
        else:
            content = f.read()
        return content

    def read_track(self, filename: str, stream: bool=False) -> GenomicTrack:
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
        return GenomicTrack.from_bedgraph(content, self._chrom_sizes)

    def __read_mask(self, filename: str , stream: bool = False) -> GenomicTrack:
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
        buffer_type = Bed6Buffer if stranded else BedBuffer
        content = self._open(filename, stream, buffer_type=buffer_type)
        return GenomicIntervals.from_intervals(content, self._chrom_sizes)

    @property
    def size(self):
        return sum(self._chrom_sizes.values())
