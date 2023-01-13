import dataclasses
from typing import List, Union, Iterable, Tuple, Dict
from .intervals import get_boolean_mask, GenomicRunLengthArray, get_pileup
from .global_offset import GlobalOffset
from ..datatypes import Interval, BedGraph
import numpy as np

GenomeIndex = Union[str, List[str], Interval, Interval.single_entry]


class MultiRLE:
    def __init__(self):
        pass

    def __getitem__(self, idx: GenomeIndex):
        pass

    @classmethod
    def from_dict(cls, d: Dict[str, GenomicRunLengthArray]) -> 'MultiRLE':
        pass

    @classmethod
    def from_stream(cls, stream: Iterable[Tuple[str, GenomicRunLengthArray]]) -> 'MultiRLE':
        pass

    @classmethod
    def from_global_pileup(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'MultiRLE':
        pass

    @classmethod
    def to_bedgraph(self) -> BedGraph:
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

    def get_pileup(self, intervals: Interval) -> GenomicRunLengthArray:
        if intervals is stream:
            groups = groupby(intervals, chromosome)
            return GroupedStream((chromosome_name, self.get_pileup(ints)) for chromosome_name, ints in groups)
        go = self._global_offset.from_local_interval(intervals)
        return get_pileup(go, self._global_size)

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
