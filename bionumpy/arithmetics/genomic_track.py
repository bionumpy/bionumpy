import numpy as np
from .global_offset import GlobalOffset
from ..computation_graph import StreamNode, Node, ComputationNode
from ..datatypes import Interval, BedGraph
from .intervals import GenomicRunLengthArray
from npstructures import RunLengthRaggedArray
from typing import List, Union, Iterable, Tuple, Dict

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

    def extract_intervals(self, intervals: Interval, stranded: bool = False):
        return self.get_intervals(intervals, stranded)

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

    def get_data(self):
        return NotImplemented        


class GenomicTrack(GenomicData):
    def sum(self) -> float:
        return NotImplemented

    def to_dict(self) -> Dict[str, GenomicRunLengthArray]:
        return NotImplemented

    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        return NotImplemented

    def to_bedgraph(self) -> 'BedGraph':
        return NotImplemented

    def to_intervals(self) -> Union[Interval, 'BedGraph']:
        pass

    @classmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'GenomicTrack':
        return GenomicTrackGlobal(global_pileup, global_offset)

    # @classmethod
    # def from_stream(cls, stream: Iterable[Tuple[str, GenomicRunLengthArray]], global_offset: GlobalOffset) -> 'GenomicTrack':
    #     return GenomicTrackStream(stream, global_offset)

    def _get_intervals_from_data(self, name, data):
        if data.dtype == bool:
            return Interval([name]*len(data.starts),
                            data.starts, data.ends)[data.values]
        else:
            return BedGraph([name]*len(data.starts),
                            data.starts, data.ends, data.values)

    @classmethod
    def from_node(cls, node):
        pass


class GenomicTrackGlobal(GenomicTrack, np.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, global_track: GenomicRunLengthArray, global_offset: GlobalOffset):
        assert isinstance(global_track, GenomicRunLengthArray)
        self._global_track = global_track
        self._global_offset = global_offset
        names = self._global_offset.names()
        starts = self._global_offset.get_size(names)
        self._genome_context = dict(zip(names, starts))

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

    # def _mask_wrapper(self, output):
    #     assert isinstance(r, GenomicRunLengthArray)
    #     if r.dtype == bool:
    #         return GenomicMaskGlobal(r, self._global_offset)
    
    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        # TODO: check that global offsets are the same
        inputs = [(i._global_track if isinstance(i, GenomicTrackGlobal) else i) for i in inputs]
        r = self._global_track.__array_ufunc__(ufunc, method, *inputs, **kwargs)
        # assert isinstance(r, GenomicRunLengthArray)
        # if r.dtype == bool:
        #    return GenomicMaskGlobal(r, self._global_offset)
        return self.__class__(r, self._global_offset)

    def get_data(self) -> Union[Interval, BedGraph]:
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
    def get_intervals(self, intervals: Union[Interval, 'GenomicIntervals'], stranded: bool = True) -> RunLengthRaggedArray:
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


class GenomicTrackNode(GenomicTrack, np.lib.mixins.NDArrayOperatorsMixin):
    def __str__(self):
        return 'GTN:' + str(self._run_length_node)

    def __init__(self, run_length_node: ComputationNode, chrom_sizes: Dict[str, int]):
        self._run_length_node = run_length_node
        self._chrom_sizes = chrom_sizes
        self._chrom_name_node = StreamNode(iter(chrom_sizes.keys()))
        self._genome_context = self._chrom_sizes

    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        args = [gtn._run_length_node if isinstance(gtn, GenomicTrackNode) else gtn for gtn in inputs]
        return self.__class__(ufunc(*args, **kwargs), self._chrom_sizes)

    def get_data(self):
        return ComputationNode(self._get_intervals_from_data, [self._chrom_name_node, self._run_length_node])

    def extract_intervals(self, intervals: Interval, stranded: bool = False):
        assert stranded is False
        return ComputationNode(lambda ra, start, stop: ra[start:stop], [self._run_length_node, intervals.start, intervals.stop])
    # GenomicRunLengthArray.extract_intervals,
    #         [self._run_length_node, intervals])
                               
# 
# class GenomicTrackStream(GenomicTrack, np.lib.mixins.NDArrayOperatorsMixin):
#     def __init__(self, track_stream: Iterable[Tuple[str, GenomicRunLengthArray]], global_offset: GlobalOffset):
#         self._track_stream = track_stream
#         self._global_offset = global_offset
# 
#     def sum(self):
#         return sum(track.sum(axis=None) for name, track in self._track_stream)
# 
#     def to_dict(self):
#         return {name: track.to_array() for name, track in self._track_stream}
# 
#     def _iter(self):
#         return self._track_stream
# 
#     def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
#         if not method == '__call__':
#             return NotImplemented
# 
#         def _args_stream(args, stream_indices):
#             args = list(args)
#             streams = tuple(args[i]._iter() for i in stream_indices)
#             for stream_args in zip(*streams):
#                 new_args = list(args)
#                 for i, (name, stream_arg) in zip(stream_indices, stream_args):
#                     new_args[i] = stream_arg
#                 yield (name, new_args)
# 
#         stream_indices = [i for i, arg in enumerate(inputs) if isinstance(arg, GenomicTrackStream)]
#         new_stream = ((name, ufunc(*new_inputs, **kwargs))
#                       for name, new_inputs in _args_stream(inputs, stream_indices))
#         return self.__class__(new_stream, self._global_offset)
# 
#         inputs = [(i._global_track if isinstance(i, GenomicTrackGlobal) else i) for i in inputs]
#         r = self._global_track.__array_ufunc__(ufunc, method, *inputs, **kwargs)
#         assert isinstance(r, GenomicRunLengthArray)
#         if r.dtype == bool:
#             return GenomicMaskGlobal(r, self._global_offset)
#         return self.__class__(r, self._global_offset)
# 
#     def to_intervals(self) -> Union[Interval, BedGraph]:
#         return NpDataclassStream(self._get_intervals_from_data(name, data)
#                                  for name, data in self._track_stream)
# 
#     def get_intervals(self, intervals: Interval, stranded: bool = True) -> RunLengthRaggedArray:
#         global_intervals = self._global_offset.from_local_interval(intervals)
#         return self._global_track[global_intervals]
# 
#     def compute(self):
#         global_track = np.concatenate([track for _, track in self._track_stream])
#         return GenomicTrackGlobal(global_track, self._global_offset)
# 
# 
# class GenomicMask(GenomicTrack):
#     @classmethod
#     def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'MultiRLE':
#         assert isinstance(global_pileup, GenomicRunLengthArray)
#         return GenomicMaskGlobal(global_pileup, global_offset)
# 
# 
# class GenomicMaskGlobal(GenomicTrackGlobal):
#     def to_intervals(self):
#         names = self._global_offset.names()
#         starts = self._global_offset.get_offset(names)
#         stops = starts+self._global_offset.get_size(names)
#         intervals_list = []
#         for name, start, stop in zip(names, starts, stops):
#             assert isinstance(self._global_track, GenomicRunLengthArray)
#             data = self._global_track[start:stop]
#             assert isinstance(data, GenomicRunLengthArray)
#             intervals_list.append(
#                 Interval([name]*len(data.starts),
#                          data.starts, data.ends)[data.values])
#         names = self._global_offset.names()
#         chrom_sizes = dict(zip(names, self._global_offset.get_size(names)))
#         return GenomicIntervals.from_intervals(
#             np.concatenate(intervals_list), chrom_sizes)
#     # return np.concatenate(intervals_list)
