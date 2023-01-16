import dataclasses
from typing import List, Union, Iterable, Tuple, Dict
from .intervals import get_boolean_mask, GenomicRunLengthArray, get_pileup
from .global_offset import GlobalOffset
from ..datatypes import Interval, BedGraph
import numpy as np

GenomeIndex = Union[str, List[str], Interval, Interval.single_entry]


class GenomicData:
    def __init__(self):
        pass

    def __getitem__(self, idx: GenomeIndex):
        pass

    def get_chromsome(self, chromosome: Union[str, List[str]]) -> 'GenomicData':
        pass

    def get_intervals(self, intervals: Interval, stranded: bool = False):
        pass

    @classmethod
    def from_dict(cls, d: Dict[str, GenomicRunLengthArray]) -> 'GenomicData':
        pass

    @classmethod
    def from_stream(cls, stream: Iterable[Tuple[str, GenomicRunLengthArray]]) -> 'MultiRLE':
        pass

    @classmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'MultiRLE':
        pass

    @classmethod
    def to_bedgraph(self) -> BedGraph:
        pass

    def to_dict(self) -> Dict[str, np.ndarray]:
        pass


class GenomicTrack(GenomicData):

    @classmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'MultiRLE':
        return GenomicTrackGlobal(global_pileup, global_offset)

    def to_dict(self):
        pass


class GenomicTrackGlobal(GenomicTrack):
    def __init__(self, global_track, global_offset):
        self._global_track = global_track
        self._global_offset = global_offset

    def to_dict(self):
        names = self._global_offset.names()
        offsets = self._global_offset.get_offset(names)
        sizes = self._global_offset.get_size(names)
        return {name: self._global_track[offset:offset+size].to_array()
                for name, offset, size in zip(names, offsets, sizes)}
    

class GenomicMask(GenomicTrack):
    pass


'''
def decorator(func):
    def new_func(intervals, *args, **kwargs):
        if intervals is stream:
            groups = groupby(intervals, 'chromosome')
            return MultiRLE.from_stream((chrom, func(i, *args, **kwargs)) for chrom, i in groups))
        return MultiRlE.from_global_offset(chrom, func)
    return new_func
'''


class Geometry:
    def __init__(self, chrom_sizes: dict):
        self._chrom_sizes = chrom_sizes
        self._global_offset = GlobalOffset(chrom_sizes)
        self._global_size = sum(chrom_sizes.values())

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
        go = self._global_offset.from_local_interval(intervals)
        return get_boolean_mask(go, self._global_size)

    def get_pileup(self, intervals: Interval) -> GenomicTrack:
        # if intervals is stream:
        #     groups = groupby(intervals, chromosome)
        #     return GroupedStream((chromosome_name, self.get_pileup(ints)) 
        #                          for chromosome_name, ints in groups)
        go = self._global_offset.from_local_interval(intervals)
        return GenomicTrack.from_global_data(
            get_pileup(go, self._global_size),
            self._global_offset)

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

    def clip(self, intervals):
        chrom_sizes = self._global_offset.get_size(intervals.chromosome)
        return dataclasses.replace(
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
        chrom_sizes = self._global_offset.get_size(intervals.chromosome)
        is_forward = intervals.strand.ravel() == "+"
        start = np.where(is_forward,
                         intervals.start,
                         np.maximum(intervals.stop-fragment_length, 0))
        stop = np.where(is_forward,
                        np.minimum(intervals.start+fragment_length, chrom_sizes),
                        intervals.stop)
                        
        return dataclasses.replace(
            intervals,
            start=start,
            stop=stop)

    def chrom_size(self, chromsome):
        return self._chrom_sizes[chromsome]

    def names(self):
        return list(self._chrom_sizes.keys())
