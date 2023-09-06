from typing import List, Union, Iterable, Tuple, Dict
from bionumpy.util.formating import table
from ..streams import groupby, NpDataclassStream
from ..streams.left_join import left_join
from ..arithmetics.intervals import get_boolean_mask, GenomicRunLengthArray, get_pileup, merge_intervals, extend_to_size
from ..datatypes import Interval, BedGraph, ChromosomeSize
from .genome_context import GenomeContext
from ..bnpdataclass import replace
from .genomic_track import GenomicArray
import numpy as np


class GeometryBase:
    def __init__(self, chrom_sizes: dict):
        self._chrom_sizes = chrom_sizes
        self._genome_context = GenomeContext.from_dict(chrom_sizes)
        # self._global_offset = GlobalOffset(chrom_sizes)
        self._global_size = self._genome_context.size

    @classmethod
    def from_chrom_sizes(cls, chrom_sizes: ChromosomeSize):
        """Create a Geometry object from a ChromosomeSizes bnpdataclass object

        Parameters
        ----------
        chrom_sizes : ChromosomeSize
        """
        return cls({str(chrom_size.name): chrom_size.size for chrom_size in chrom_sizes})

    def chrom_size(self, chromsome: str) -> int:
        """Return the size of the given chromosome

        Parameters
        ----------
        chromsome : str

        Returns
        -------
        int

        """
        return self._chrom_sizes[chromsome]

    def names(self) -> List[str]:
        """List the chromosomes in order

        Returns
        -------
        List[str]

        """

        return list(self._chrom_sizes.keys())

    def size(self) -> int:
        """Return the total size of the genome

        Returns
        -------
        int
        """
        return self._genome_context.size


class Geometry(GeometryBase):
    def jaccard(self, intervals_a: Interval, intervals_b: Interval) -> float:
        """Calculate the Jaccard similarity score between two sets of intervals

        Parameters
        ----------
        intervals_a : Interval
        intervals_b : Interval

        Returns
        -------
        float
            Jaccard similarity score

        """
        a = self.get_global_mask(intervals_a)
        b = self.get_global_mask(intervals_b)
        intersect = (a & b).sum()
        union = a.sum() + b.sum() - intersect
        assert union >= 0
        return intersect / union

    def get_global_mask(self, intervals: Interval) -> GenomicRunLengthArray:
        """Create a mask of any area covered by any interval along the entire genome

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        GenomicRunLengthArray
            Runlength encoded boolean mask

        """
        if isinstance(intervals, GenomicRunLengthArray):
            return intervals
        go = self._genome_context.global_offset.from_local_interval(intervals)
        return get_boolean_mask(go, self._genome_context.size)

    def get_mask(self, intervals: Interval) -> 'GenomicMask':
        """Create a GenomeMask of all areas covered by at least one interval

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        GenomicMask
        """
        return GenomicArray.from_global_data(self.get_global_mask(intervals), self._genome_context)

    def get_pileup(self, intervals: Interval) -> 'GenomicArray':
        """Create a GenomicTrack of how many intervals covers each position in the genome

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        GenomicArray
        """
        go = self._genome_context.global_offset.from_local_interval(intervals)
        return GenomicArray.from_global_data(
            get_pileup(go, self._global_size),
            self._genome_context)

    def jaccard_all_vs_all(self, intervals_list: List[Interval]) -> np.ndarray:
        """Calculate all pairwise Jaccard similarity scores for a list of interval sets

        Parameters
        ----------
        intervals_list : List[Interval]

        Returns
        -------
        np.ndarray
            Matrix containing all pairwise Jaccard similarity scores

        """
        masks = [self.get_global_mask(intervals) for intervals in intervals_list]
        matrix = np.zeros((len(intervals_list), len(intervals_list)))
        for i, a in enumerate(masks):
            for j, b in enumerate(masks[i+1:], 1):
                result = self.jaccard(a, b)
                matrix[i, i+j] = result
                matrix[i+j, i] = result

        return matrix

    def clip(self, intervals: Interval) -> Interval:
        """Clip intervals so that all intervals are contained in their corresponding chromosome

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        Interval
        """
        chrom_sizes = self._genome_context.global_offset.get_size(intervals.chromosome)
        return replace(
            intervals,
            start=np.maximum(0, intervals.start),
            stop=np.minimum(chrom_sizes, intervals.stop))

    def extend_to_size(self, intervals: Interval, fragment_length: int) -> Interval:
        """Extend/shrink intervals to match the given length. Stranded

        For intervals on the + strand, keep the start coordinate and
        adjust the stop coordinate. For - intervals keep the stop
        coordinate and adjust the start

        Parameters
        ----------
        intervals : Interval
        fragment_length : int

        Returns
        -------
        Interval
        """
        chrom_sizes = self._genome_context.global_offset.get_size(intervals.chromosome)
        return extend_to_size(intervals, fragment_length, chrom_sizes)

    def get_track(self, bedgraph: BedGraph) -> GenomicArray:
        """Create a genomic track from a bedgraph

        Parameters
        ----------
        bedgraph : BedGraph

        Returns
        -------
        GenomicArray
        """
        gi = self._genome_context.global_offset.from_local_interval(bedgraph)
        rle = GenomicRunLengthArray.from_bedgraph(gi)
        return GenomicArray.from_global_data(rle, self._genome_context)

    def merge_intervals(self, intervals: Interval, distance: int = 0) -> Interval:
        """Merge all intervals that either overlaps or lie within a given distance

        Parameters
        ----------
        intervals : Interval
        distance : int

        Returns
        -------
        Interval

        """
        global_intervals = self._genome_context._global_offset.from_local_interval(intervals)
        global_merged = merge_intervals(global_intervals, distance)
        return self._genome_context._global_offset.to_local_interval(global_merged)

    def sort(self, intervals: Interval) -> Interval:
        """Sort a set of intervals according to the chormosome order

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        Interval

        """
        global_intervals = self._genome_context.global_offset.from_local_interval(intervals)
        return self._genome_context.global_offset.to_local_interval(global_intervals.sort_by('start'))

    def __repr__(self):
        return f"{self.__class__.__name__}(" + repr(self._chrom_sizes) + ")"

    def __str__(self):
        return table(zip(self._chrom_sizes.keys(), self._chrom_sizes.values()), headers=["Chromosome", "Size"])


class StreamedGeometry(GeometryBase):
    def get_track(self, bedgraph: Iterable[BedGraph]) -> GenomicArray:
        grouped = groupby(bedgraph, 'chromosome')
        track_stream = ((name, GenomicRunLengthArray.from_bedgraph(b))
                        for name, b in grouped)
        return GenomicArray.from_stream(track_stream, self._genome_context)

    def get_pileup(self, intervals: Iterable[Interval]) -> GenomicArray:
        grouped = groupby(intervals, 'chromosome')
        pileups = ((name, get_pileup(intervals, size)) if intervals is not None else GenomicRunLengthArray(np.array([0, size]), [0]) for name, size, intervals in left_join(self._chrom_sizes.items(), grouped))
        return GenomicArray.from_stream(pileups, self._genome_context)

    def extend_to_size(self, intervals: Iterable[Interval], fragment_length: int) -> Iterable[Interval]:
        """Extend/shrink intervals to match the given length. Stranded

        For intervals on the + strand, keep the start coordinate and
        adjust the stop coordinate. For - intervals keep the stop
        coordinate and adjust the start

        Parameters
        ----------
        intervals : Interval
        fragment_length : int

        Returns
        -------
        Interval
        """
        if isinstance(intervals, GenomicIntervalsStreamed):
            return Geometry(self._chrom_sizes).extend_to_size(intervals, fragment_length)

        return NpDataclassStream(Geometry(self._chrom_sizes).extend_to_size(i, fragment_length)
                                 for i in intervals)

    def clip(self, intervals: Interval) -> Interval:
        """Clip intervals so that all intervals are contained in their corresponding chromosome

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        Interval
        """
        return NpDataclassStream(Geometry(self._chrom_sizes).clip(i)
                                 for i in intervals)

    def merge_intervals(self, intervals: Iterable[Interval], distance: int = 0) -> Iterable[Interval]:
        """Merge all intervals that either overlaps or lie within a given distance

        Parameters
        ----------
        intervals : Interval
        distance : int

        Returns
        -------
        Interval

        """
        return NpDataclassStream(merge_intervals(i, distance) for i in intervals)
