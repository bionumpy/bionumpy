import numpy as np
from abc import abstractclassmethod, abstractmethod, abstractproperty, ABC
from .global_offset import GlobalOffset
from ..computation_graph import StreamNode, Node, ComputationNode
from ..datatypes import Interval, BedGraph
from ..arithmetics.intervals import GenomicRunLengthArray
from ..bnpdataclass import BNPDataClass
from ..streams import groupby
from npstructures import RunLengthRaggedArray
from typing import List, Union, Iterable, Tuple, Dict

GenomeIndex = Union[str, List[str], Interval, Interval.single_entry]


class GenomicData:

    def __getitem__(self, idx: GenomeIndex):
        if isinstance(idx, str):
            return self.extract_chromsome(idx)
        if isinstance(idx, Interval) or (hasattr(idx, 'start') and hasattr(idx, 'stop') and hasattr(idx, 'chromosome')):
            return self.extract_intervals(idx)
        if isinstance(idx, list):
            if len(idx) == 0:
                return self.empty()
            if isinstance(idx[0], str):
                return self.extract_intervals(idx)
        assert False

    @abstractmethod
    def extract_chromsome(self, chromosome: Union[str, List[str]]) -> 'GenomicData':
        return NotImplemented

    @abstractmethod
    def extract_intervals(self, intervals: Interval, stranded: bool = False) -> RunLengthRaggedArray:
        """Get the data within the (stranded) intervals

        Parameters
        ----------
        intervals : Interval
            Set of intervals
        stranded : bool
            Wheter to reverse intervals on - strand

        Returns
        -------
        RunLengthRaggedArray
            Data for all intervals
        """
        return NotImplemented

    @abstractclassmethod
    def from_dict(cls, d: Dict[str, GenomicRunLengthArray]) -> 'GenomicData':
        """Create genomic data from a dict of data with chromosomes as keys

        Parameters
        ----------
        d : Dict[str, GenomicRunLengthArray]

        Returns
        -------
        'GenomicData'
        """
        
        return NotImplemented

    @abstractclassmethod
    def from_stream(cls, stream: Iterable[Tuple[str, GenomicRunLengthArray]], chrom_sizes: dict) -> 'GenomicData':
        return NotImplemented

    @abstractclassmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'GenomicData':
        return NotImplempented

    @abstractmethod
    def to_dict(self) -> Dict[str, np.ndarray]:
        return NotImplemented

    @abstractmethod
    def get_data(self):
        return NotImplemented


class GenomicTrack(GenomicData):

    @abstractmethod
    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        return NotImplemented

    def __array_function__(self, func: callable, types: List, args: List, kwargs: Dict):
        return NotImplemented

    @abstractmethod
    def sum(self, axis=None) -> float:
        return NotImplemented

    def to_bedgraph(self) -> 'BedGraph':
        return NotImplemented

    @classmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'GenomicTrack':
        return GenomicTrackGlobal(global_pileup, global_offset)

    def _get_intervals_from_data(self, name, data):
        if data.dtype == bool:
            return Interval([name]*len(data.starts),
                            data.starts, data.ends)[data.values]
        else:
            return BedGraph([name]*len(data.starts),
                            data.starts, data.ends, data.values)

    @classmethod
    def from_bedgraph(cls, bedgraph: BedGraph, chrom_sizes: Dict[str, int]):
        if isinstance(bedgraph, BedGraph):
            go = GlobalOffset(chrom_sizes)
            gi = go.from_local_interval(bedgraph)
            rle = GenomicRunLengthArray.from_bedgraph(gi, go.total_size())
            return cls.from_global_data(rle, go)

        filled = fill_grouped(groupby(bedgraph, 'chromosome'), chrom_sizes.keys(), BedGraph)
        interval_stream = StreamNode(filled)
        return GenomicTrackNode(ComputationNode(GenomicRunLengthArray.from_bedgraph,
                                                [interval_stream, StreamNode(iter(chrom_sizes.values()))]),
                                chrom_sizes)


def fill_grouped(grouped, real_order, dataclass):
    real_order = iter(real_order)
    next_real = next(real_order, None)
    for name, group in grouped:
        assert next_real is not None
        while name != next_real:
            yield dataclass.empty()
            next_real = next(real_order, None)
        yield group
        next_real = next(real_order, None)
    # next_real = next(real_order, None)
    while next_real is not None:
        yield next_real
        next_real = next(real_order, None)


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

    def extract_chromsome(self, chromosome):
        assert isinstance(chromosome, str)
        offset = self._global_offset.get_offset([chromosome])[0]
        return self._global_track[offset:offset + self._global_offset.get_size([chromosome])[0]]

    def __str__(self) -> str:
        return str(self._global_track)

    def to_dict(self) -> Dict[str, GenomicRunLengthArray]:
        names = self._global_offset.names()
        offsets = self._global_offset.get_offset(names)
        sizes = self._global_offset.get_size(names)
        return {name: self._global_track[offset:offset+size].to_array()
                for name, offset, size in zip(names, offsets, sizes)}

    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        inputs = [(i._global_track if isinstance(i, GenomicTrackGlobal) else i) for i in inputs]
        r = self._global_track.__array_ufunc__(ufunc, method, *inputs, **kwargs)
        return self.__class__(r, self._global_offset)

    def __array_function__(self, func: callable, types: List, args: List, kwargs: Dict):
        """Handles any numpy array functions called on a raggedarray

        Parameters
        ----------
        func : callable
        types : List
        args : List
        kwargs : Dict
        """
        args = [(i._global_track if isinstance(i, GenomicTrackGlobal) else i) for i in args]
        if func == np.histogram:
            return np.histogram(*args, **kwargs)

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

    def extract_intervals(self, intervals: Union[Interval, 'GenomicIntervals'], stranded: bool = False) -> RunLengthRaggedArray:
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

    def extract_intervals(self, intervals: Interval, stranded: bool = False) -> RunLengthRaggedArray:
        assert stranded is False
        return ComputationNode(lambda ra, start, stop: ra[start:stop], [self._run_length_node, intervals.start, intervals.stop])

    def __array_function__(self, func: callable, types: List, args: List, kwargs: Dict):
        """Handles any numpy array functions called on a raggedarray

        Parameters
        ----------
        func : callable
        types : List
        args : List
        kwargs : Dict
        """
        if func == np.histogram:
            return np.histogram(args[0]._run_length_node, *args[1:], **kwargs)
        return NotImplemented
