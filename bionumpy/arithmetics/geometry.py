import dataclasses
from typing import List, Union, Iterable, Tuple, Dict
from bionumpy.util.formating import table
from ..streams import groupby, NpDataclassStream
from ..streams.left_join import left_join
from .intervals import get_boolean_mask, GenomicRunLengthArray, get_pileup, merge_intervals
from .global_offset import GlobalOffset
from ..datatypes import Interval, BedGraph, ChromosomeSize
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
        return NotImplemented

    @classmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'MultiRLE':
        return NotImplemented

    def to_dict(self) -> Dict[str, np.ndarray]:
        return NotImplemented

    def compute(self):
        return NotImplemented


class GenomicTrack(GenomicData):
    def sum(self) -> float:
        return NotImplemented

    def to_dict(self) -> Dict[str, GenomicRunLengthArray]:
        return NotImplemented

    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        return NotImplemented

    def to_bedgraph(self) -> BedGraph:
        return NotImplemented

    def to_intervals(self) -> Union[Interval, BedGraph]:
        pass

    @classmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'GenomicTrack':
        return GenomicTrackGlobal(global_pileup, global_offset)

    @classmethod
    def from_stream(cls, stream: Iterable[Tuple[str, GenomicRunLengthArray]], global_offset: GlobalOffset) -> 'GenomicTrack':
        return GenomicTrackStream(stream, global_offset)

    def _get_intervals_from_data(name, data):
        if data.dtype == bool:
            return Interval([name]*len(data.starts),
                            data.starts, data.ends)[data.values]
        else:
            return BedGraph([name]*len(data.starts),
                            data.starts, data.ends, data.values)


class GenomicTrackGlobal(GenomicTrack, np.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, global_track: GenomicRunLengthArray, global_offset: GlobalOffset):
        assert isinstance(global_track, GenomicRunLengthArray)
        self._global_track = global_track
        self._global_offset = global_offset

    def sum(self) -> float:
        return self._global_track.sum(axis=None)

    def __str__(self) -> str:
        return str(self._global_track)

    def to_dict(self) -> Dict[str, GenomicRunLengthArray]:
        names = self._global_offset.names()
        offsets = self._global_offset.get_offset(names)
        sizes = self._global_offset.get_size(names)
        return {name: self._global_track[offset:offset+size].to_array()
                for name, offset, size in zip(names, offsets, sizes)}

    def _mask_wrapper(self, output):
        assert isinstance(r, GenomicRunLengthArray)
        if r.dtype == bool:
            return GenomicMaskGlobal(r, self._global_offset)
    
    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        # TODO: check that global offsets are the same
        inputs = [(i._global_track if isinstance(i, GenomicTrackGlobal) else i) for i in inputs]
        r = self._global_track.__array_ufunc__(ufunc, method, *inputs, **kwargs)
        # assert isinstance(r, GenomicRunLengthArray)
        if r.dtype == bool:
            return GenomicMaskGlobal(r, self._global_offset)
        return self.__class__(r, self._global_offset)

    def to_intervals(self) -> Union[Interval, BedGraph]:
        names = self._global_offset.names()
        starts = self._global_offset.get_offset(names)
        stops = starts+self._global_offset.get_size(names)
        intervals_list = []
        for name, start, stop in zip(names, starts, stops):
            assert isinstance(self._global_track, GenomicRunLengthArray)
            data = self._global_track[start:stop]
            assert isinstance(data, GenomicRunLengthArray)
            intervals_list.append(self._get_intervals_from_data(name, data))
        return np.concatenate(intervals_list)

    def to_bedgraph(self) -> BedGraph:
        """Create a BedGraph object corresponding to this track

        Returns
        -------
        BedGraph
        """
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

    # extract
    def get_intervals(self, intervals: Interval, stranded: bool = True) -> RunLengthRaggedArray:
        """Extract the data contained in a set of intervals

        Parameters
        ----------
        intervals : Interval
        stranded : bool

        Returns
        -------
        RunLengthRaggedArray
        """
        
        global_intervals = self._global_offset.from_local_interval(intervals)
        return self._global_track[global_intervals]

    def compute(self) -> GenomicTrack:
        return self


class GenomicTrackStream(GenomicTrack, np.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, track_stream: Iterable[Tuple[str, GenomicRunLengthArray]], global_offset: GlobalOffset):
        self._track_stream = track_stream
        self._global_offset = global_offset

    def sum(self):
        return sum(track.sum(axis=None) for name, track in self._track_stream)

    def to_dict(self):
        return {name: track.to_array() for name, track in self._track_stream}

    def _iter(self):
        return self._track_stream

    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        if not method == '__call__':
            return NotImplemented

        def _args_stream(args, stream_indices):
            args = list(args)
            streams = tuple(args[i]._iter() for i in stream_indices)
            for stream_args in zip(*streams):
                new_args = list(args)
                for i, (name, stream_arg) in zip(stream_indices, stream_args):
                    new_args[i] = stream_arg
                yield (name, new_args)

        stream_indices = [i for i, arg in enumerate(inputs) if isinstance(arg, GenomicTrackStream)]
        new_stream = ((name, ufunc(*new_inputs, **kwargs))
                      for name, new_inputs in _args_stream(inputs, stream_indices))
        return self.__class__(new_stream, self._global_offset)
        inputs = [(i._global_track if isinstance(i, GenomicTrackGlobal) else i) for i in inputs]
        r = self._global_track.__array_ufunc__(ufunc, method, *inputs, **kwargs)
        assert isinstance(r, GenomicRunLengthArray)
        if r.dtype == bool:
            return GenomicMaskGlobal(r, self._global_offset)
        return self.__class__(r, self._global_offset)

    def to_intervals(self) -> Union[Interval, BedGraph]:
        return NpDataclassStream(self._get_intervals_from_data(name, data)
                                 for name, data in self._track_stream)

    def get_intervals(self, intervals: Interval, stranded: bool = True) -> RunLengthRaggedArray:
        global_intervals = self._global_offset.from_local_interval(intervals)
        return self._global_track[global_intervals]

    def compute(self):
        global_track = np.concatenate([track for _, track in self._track_stream])
        return GenomicTrackGlobal(global_track, self._global_offset)


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


class GeometryBase:
    def __init__(self, chrom_sizes: dict):
        self._chrom_sizes = chrom_sizes
        self._global_offset = GlobalOffset(chrom_sizes)
        self._global_size = sum(chrom_sizes.values())

    @classmethod
    def from_chrom_sizes(cls, chrom_sizes: ChromosomeSize):
        """Create a Geometry object from a ChromosomeSizes bnpdataclass object

        Parameters
        ----------
        chrom_sizes : ChromosomeSize
        """
        return cls({chrom_size.name.to_string():chrom_size.size
                    for chrom_size in chrom_sizes})

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
        return self._global_size



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
        go = self._global_offset.from_local_interval(intervals)
        return get_boolean_mask(go, self._global_size)

    def get_mask(self, intervals: Interval) -> GenomicMask:
        """Create a GenomeMask of all areas covered by at least one interval

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        GenomicMask
        """
        return GenomicMask.from_global_data(self.get_global_mask(intervals), self._global_offset)

    def get_pileup(self, intervals: Interval) -> GenomicTrack:
        """Create a GenomicTrack of how many intervals covers each position in the genome

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        GenomicTrack
        """
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

    def clip(self, intervals: Interval) -> Interval:
        """Clip intervals so that all intervals are contained in their corresponding chromosome

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        Interval
        """
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
        """Create a genomic track from a bedgraph

        Parameters
        ----------
        bedgraph : BedGraph

        Returns
        -------
        GenomicTrack
        """
        gi = self._global_offset.from_local_interval(bedgraph)
        rle = GenomicRunLengthArray.from_bedgraph(gi)
        return GenomicTrack.from_global_data(rle, self._global_offset)

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
        global_intervals = self._global_offset.from_local_interval(intervals)
        global_merged = merge_intervals(global_intervals, distance)
        return self._global_offset.to_local_interval(global_merged)

    def sort(self, intervals: Interval) -> Interval:
        """Sort a set of intervals according to the chormosome order

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        Interval

        """
        global_intervals = self._global_offset.from_local_interval(intervals)
        return self._global_offset.to_local_interval(global_intervals.sort_by('start'))

    def __repr__(self):
        return f"{self.__class__.__name__}(" + repr(self._chrom_sizes) + ")"

    def __str__(self):
        return table(zip(self._chrom_sizes.keys(), self._chrom_sizes.values()), headers=["Chromosome", "Size"])



class StreamedGeometry(GeometryBase):
    def get_track(self, bedgraph: Iterable[BedGraph]) -> GenomicTrack:
        grouped = groupby(bedgraph, 'chromosome')
        track_stream = ((name, GenomicRunLengthArray.from_bedgraph(b))
                        for name, b in grouped)
        return GenomicTrack.from_stream(track_stream, self._global_offset)

    def get_pileup(self, intervals: Iterable[Interval]) -> GenomicTrack:
        grouped = groupby(intervals, 'chromosome')
        pileups = ((name, get_pileup(intervals, size)) if intervals is not None else GenomicRunLengthArray(np.array([0, size]), [0]) for name, size, intervals in left_join(self._chrom_sizes.items(), grouped))
        return GenomicTrack.from_stream(pileups, self._global_offset)
            
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
        print('Here')
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
        return self._global_size
