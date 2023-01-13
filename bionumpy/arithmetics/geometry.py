import dataclasses
from .intervals import get_boolean_mask, GenomicRunLengthArray
from .global_offset import GlobalOffset
from ..datatypes import Interval
import numpy as np


class Geometry:
    def __init__(self, chrom_sizes: dict):
        self._chrom_sizes = chrom_sizes
        self._global_offset = GlobalOffset(chrom_sizes)
        self._global_size = sum(chrom_sizes.values())

    def jaccard(self, intervals_a: Interval, intervals_b: Interval) -> float:
        a = self.get_global_mask(intervals_a)
        b = self.get_global_mask(intervals_b)
        intersect = (a & b).sum()
        union = a.sum() + b.sum() - intersect
        assert union >= 0
        return intersect / union

    def get_global_mask(self, intervals: Interval):
        if isinstance(intervals, GenomicRunLengthArray):
            return intervals
        go = self._global_offset.from_local_interval(intervals)
        return get_boolean_mask(go, self._global_size)

    def jaccard_all_vs_all(self, intervals_list):
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
