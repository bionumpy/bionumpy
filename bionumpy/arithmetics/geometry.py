import dataclasses
from typing import List, Union, Iterable, Tuple, Dict
from ..streams import groupby
from .intervals import get_boolean_mask, GenomicRunLengthArray, get_pileup, merge_intervals
from .global_offset import GlobalOffset
from ..datatypes import Interval, BedGraph
from npstructures import RunLengthRaggedArray
import numpy as np

GenomeIndex = Union[str, List[str], Interval, Interval.single_entry]


class GenomicData:
    def __init__(self):
        pass

    def __getitem__(self, idx: GenomeIndex):
        return NotImplemented

    def get_chromsome(self, chromosome: Union[str, List[str]]) -> 'GenomicData':
        return NotImplemented

    def get_intervals(self, intervals: Interval, stranded: bool = False):
        return NotImplemented

    @classmethod
    def from_dict(cls, d: Dict[str, GenomicRunLengthArray]) -> 'GenomicData':
        return NotImplemented

    @classmethod
    def from_stream(cls, stream: Iterable[Tuple[str, GenomicRunLengthArray]], chrom_sizes: dict) -> 'MultiRLE':
        pass

    @classmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'MultiRLE':
        pass

    def to_dict(self) -> Dict[str, np.ndarray]:
        pass


class GenomicTrack(GenomicData):
    def sum(self):
        return NotImplemented

    def to_dict(self):
        return NotImplemented

    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        return NotImplemented

    def to_bedgraph(self) -> BedGraph:
        return NotImplemented

    @classmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'MultiRLE':
        return GenomicTrackGlobal(global_pileup, global_offset)

    @classmethod
    def from_stream(cls, stream: Iterable[Tuple[str, GenomicRunLengthArray]], chrom_sizes: dict) -> 'MultiRLE':
        return GenomicTrackStream(stream)

    def to_dict(self):
        pass


class GenomicTrackGlobal(GenomicTrack, np.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, global_track, global_offset):
        assert isinstance(global_track, GenomicRunLengthArray)
        self._global_track = global_track
        self._global_offset = global_offset

    def sum(self):
        return self._global_track.sum(axis=None)

    def __str__(self):
        return str(self._global_track)

    def to_dict(self):
        names = self._global_offset.names()
        offsets = self._global_offset.get_offset(names)
        sizes = self._global_offset.get_size(names)
        return {name: self._global_track[offset:offset+size].to_array()
                for name, offset, size in zip(names, offsets, sizes)}

    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        # TODO: check that global offsets are the same
        inputs = [(i._global_track if isinstance(i, GenomicTrackGlobal) else i) for i in inputs]
        r = self._global_track.__array_ufunc__(ufunc, method, *inputs, **kwargs)
        assert isinstance(r, GenomicRunLengthArray)
        if r.dtype == bool:
            return GenomicMaskGlobal(r, self._global_offset)
        return self.__class__(r, self._global_offset)

    def to_bedgraph(self) -> BedGraph:
        names = self._global_offset.names()
        starts = self._global_offset.get_offset(names)
        stops = starts+self._global_offset.get_size(names)
        intervals_list = []
        for name, start, stop in zip(names, starts, stops):
            assert isinstance(self._global_track, GenomicRunLengthArray)
            data = self._global_track[start:stop]
            assert isinstance(data, GenomicRunLengthArray)
            intervals_list.append(
                BedGraph([name]*len(data.starts),
                         data.starts, data.ends, data.values))
        return np.concatenate(intervals_list)

    def get_intervals(self, intervals: Interval, stranded: bool = True) -> RunLengthRaggedArray:
        global_intervals = self._global_offset.from_local_interval(intervals)
        return self._global_track[global_intervals]


class GenomicTrackStream(GenomicTrack):
    def __init__(self, track_stream):
        self._track_stream = track_stream

    def sum(self):
        return sum(track.sum(axis=None) for name, track in self._track_stream)

    def to_dict(self):
        return dict(self._track_stream)

    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        # TODO: check that global offsets are the same
        inputs = [(i._global_track if isinstance(i, GenomicTrackGlobal) else i) for i in inputs]
        r = self._global_track.__array_ufunc__(ufunc, method, *inputs, **kwargs)
        assert isinstance(r, GenomicRunLengthArray)
        if r.dtype == bool:
            return GenomicMaskGlobal(r, self._global_offset)
        return self.__class__(r, self._global_offset)

    def to_bedgraph(self) -> BedGraph:
        names = self._global_offset.names()
        starts = self._global_offset.get_offset(names)
        stops = starts+self._global_offset.get_size(names)
        intervals_list = []
        for name, start, stop in zip(names, starts, stops):
            assert isinstance(self._global_track, GenomicRunLengthArray)
            data = self._global_track[start:stop]
            assert isinstance(data, GenomicRunLengthArray)
            intervals_list.append(
                BedGraph([name]*len(data.starts),
                         data.starts, data.ends, data.values))
        return np.concatenate(intervals_list)

    def get_intervals(self, intervals: Interval, stranded: bool = True) -> RunLengthRaggedArray:
        global_intervals = self._global_offset.from_local_interval(intervals)
        return self._global_track[global_intervals]


class GenomicMask(GenomicTrack):
    @classmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'MultiRLE':
        assert isinstance(global_pileup, GenomicRunLengthArray)
        return GenomicMaskGlobal(global_pileup, global_offset)


class GenomicMaskGlobal(GenomicTrackGlobal):
    def to_intervals(self):
        names = self._global_offset.names()
        starts = self._global_offset.get_offset(names)
        stops = starts+self._global_offset.get_size(names)
        intervals_list = []
        for name, start, stop in zip(names, starts, stops):
            assert isinstance(self._global_track, GenomicRunLengthArray)
            data = self._global_track[start:stop]
            assert isinstance(data, GenomicRunLengthArray)
            intervals_list.append(
                Interval([name]*len(data.starts),
                         data.starts, data.ends)[data.values])
        return np.concatenate(intervals_list)


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

    def get_mask(self, intervals: Interval) -> GenomicMask:
        return GenomicMask.from_global_data(self.get_global_mask(intervals), self._global_offset)

    def get_pileup(self, intervals: Interval) -> GenomicTrack:
        #if isinstance(intervals, BnpStream):
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

    def get_track(self, bedgraph: BedGraph) -> GenomicTrack:
        gi = self._global_offset.from_local_interval(bedgraph)
        rle = GenomicRunLengthArray.from_bedgraph(gi)
        return GenomicTrack.from_global_data(rle, self._global_offset)

    def merge_intervals(self, intervals: Interval, distance: int = 0) -> Interval:
        global_intervals = self._global_offset.from_local_interval(intervals)
        global_merged = merge_intervals(global_intervals, distance)
        return self._global_offset.to_local_interval(global_merged)

    def sort(self, intervals: Interval) -> Interval:
        global_intervals = self._global_offset.from_local_interval(intervals)
        return self._global_offset.to_local_interval(global_intervals.sort_by('start'))

    def chrom_size(self, chromsome):
        return self._chrom_sizes[chromsome]

    def names(self):
        return list(self._chrom_sizes.keys())

    def size(self):
        return self._global_size

class StreamedGeometry(Geometry):
    def __init__(self, chrom_sizes: dict):
        self._chrom_sizes = chrom_sizes
        self._global_offset = GlobalOffset(chrom_sizes)
        self._global_size = sum(chrom_sizes.values())

    def get_track(self, bedgraph: BedGraph) -> GenomicTrack:
        grouped = groupby(bedgraph, 'chromosome')
        track_stream = ((name, GenomicRunLengthArray.from_bedgraph(b))
                        for name, b in grouped)
        return GenomicTrack.from_stream(track_stream, self._chrom_sizes)
