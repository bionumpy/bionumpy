import numpy as np
from abc import abstractclassmethod, abstractmethod, abstractproperty, ABC
from .global_offset import GlobalOffset
from .genomic_data import GenomicData
from ..computation_graph import StreamNode, Node, ComputationNode
from ..datatypes import Interval, BedGraph
from ..arithmetics.intervals import GenomicRunLengthArray
from ..bnpdataclass import BNPDataClass
from ..streams import groupby
from npstructures import RunLengthRaggedArray
from typing import List, Union, Iterable, Tuple, Dict


class GenomicArray(GenomicData):

    def __getitem__(self, idx):
        if isinstance(idx, GenomicArray):
            return self._index_boolean(idx)
        return super().__getitem__(idx)

    @abstractmethod
    def _index_boolean(self, idx):
        return NotImplemented

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
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, global_offset: GlobalOffset) -> 'GenomicArray':
        return GenomicArrayGlobal(global_pileup, global_offset)

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
        return GenomicArrayNode(ComputationNode(GenomicRunLengthArray.from_bedgraph,
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


class GenomicArrayGlobal(GenomicArray, np.lib.mixins.NDArrayOperatorsMixin):
    def __init__(self, global_track: GenomicRunLengthArray, global_offset: GlobalOffset):
        assert isinstance(global_track, GenomicRunLengthArray)
        self._global_track = global_track
        self._global_offset = global_offset
        names = self._global_offset.names()
        starts = self._global_offset.get_size(names)
        self._genome_context = dict(zip(names, starts))

    def _index_boolean(self, idx):
        assert idx.dtype == bool
        assert isinstance(idx, GenomicArrayGlobal)
        return self._global_track[idx._global_track]

    def sum(self) -> float:
        return self._global_track.sum(axis=None)

    def extract_chromsome(self, chromosome):
        assert isinstance(chromosome, str)
        offset = self._global_offset.get_offset([chromosome])[0]
        return self._global_track[offset:offset + self._global_offset.get_size([chromosome])[0]]

    def __repr__(self) -> str:
        lines = []
        for name, _ in zip(self._genome_context, range(10)):
            lines.append(f'{name}: {self[name]}')
        if len(self._genome_context) > 10:
            lines.extend('...')
        return '\n'.join(lines)

    def to_dict(self) -> Dict[str, GenomicRunLengthArray]:
        names = self._global_offset.names()
        offsets = self._global_offset.get_offset(names)
        sizes = self._global_offset.get_size(names)
        return {name: self._global_track[offset:offset+size].to_array()
                for name, offset, size in zip(names, offsets, sizes)}

    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        inputs = [(i._global_track if isinstance(i, GenomicArrayGlobal) else i) for i in inputs]
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
        args = [(i._global_track if isinstance(i, GenomicArrayGlobal) else i) for i in args]
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
        rle = self._global_track[global_intervals]
        if not stranded:
            return rle
        r = rle[:, ::-1]
        return np.where((intervals.strand.ravel() == '+')[:, np.newaxis],
                        rle, r)


class GenomicArrayNode(GenomicArray, np.lib.mixins.NDArrayOperatorsMixin):
    def __str__(self):
        return 'GTN:' + str(self._run_length_node)

    def __init__(self, run_length_node: ComputationNode, chrom_sizes: Dict[str, int]):
        self._run_length_node = run_length_node
        self._chrom_sizes = chrom_sizes
        self._chrom_name_node = StreamNode(iter(chrom_sizes.keys()))
        self._genome_context = self._chrom_sizes

    def __array_ufunc__(self, ufunc: callable, method: str, *inputs, **kwargs):
        args = [gtn._run_length_node if isinstance(gtn, GenomicArrayNode) else gtn for gtn in inputs]
        return self.__class__(ufunc(*args, **kwargs), self._chrom_sizes)

    def get_data(self):
        return ComputationNode(self._get_intervals_from_data, [self._chrom_name_node, self._run_length_node])

    def _index_boolean(self, idx):
        assert idx.dtype == bool
        assert isinstance(idx, GenomicArrayNode)
        return ComputationNode(lambda a, i: a[i], [self._run_length_node, idx])

    def _extract_full_intervals(self, intervals, stranded):
        return NotImplemented
        #stream_intervals = intervals.get_sorted_stream()
        #subset = self.extract_intervals(stream_intervals, stranded)
        ## return ComputationNode(subset.__getitem__(
        #return ComputationNode(lambda ra, start, stop: ra[start:stop], [self._run_length_node, intervals.start, intervals.stop])
        #
        #sorted_intervals = intervals[argsort]
        # grouped = StreamNode((g[1] for g in groupby(sorted_intervals, 'chromosome')))
        

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
