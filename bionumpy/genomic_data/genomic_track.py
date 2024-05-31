import numpy as np
from abc import abstractmethod
from .genomic_data import GenomicData
from .genome_context_base import GenomeContextBase
from ..computation_graph import StreamNode, ComputationNode
from ..datatypes import Interval, BedGraph
from ..arithmetics.intervals import GenomicRunLengthArray
from npstructures import RunLengthRaggedArray
from typing import List, Union, Iterable, Tuple, Dict, Any, Callable


class GenomicArray(GenomicData):
    '''
    Class for representing data on a genome.
    '''

    @property
    def genome_context(self) -> GenomeContextBase:
        return self._genome_context

    @abstractmethod
    def __array_ufunc__(self, ufunc: Callable, method: str, *inputs, **kwargs) -> 'GenomicArray':
        return NotImplemented

    @abstractmethod
    def __array_function__(self, func: Callable, types: List, args: List, kwargs: Dict) -> Any:
        return NotImplemented

    @abstractmethod
    def sum(self, axis: int = None) -> float:
        return NotImplemented

    def to_bedgraph(self) -> BedGraph:
        return NotImplemented

    @classmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray,
                         genome_context: GenomeContextBase) -> 'GenomicArray':
        """Create the genomic array from data represented on a flat/concatenated genome.

        Parameters
        ----------
        global_pileup : GenomicRunLengthArray
            Array over the concatenated genome
        genome_context : GenomeContextBase
            The genome context

        Returns
        -------
        'GenomicArray'
        """

        return GenomicArrayGlobal(global_pileup, genome_context)

    @classmethod
    def from_bedgraph(cls, bedgraph: BedGraph, genome_context: GenomeContextBase) -> 'GenomicData':
        """Create a genomic array from a bedgraph and genomic context

        Parameters
        ----------
        bedgraph : BedGraph
            Bedgraph with data along the genome
        genome_context : GenomeContextBase
            The genome context

        Returns
        -------
        'GenomicData'
        """

        if isinstance(bedgraph, BedGraph):
            go = genome_context.global_offset
            gi = go.from_local_interval(bedgraph)
            rle = GenomicRunLengthArray.from_bedgraph(gi, go.total_size())
            return cls.from_global_data(rle, genome_context)
        filled = genome_context.iter_chromosomes(bedgraph, BedGraph)
        # filled = fill_grouped(groupby(bedgraph, 'chromosome'), chrom_sizes.keys(), BedGraph)
        interval_stream = StreamNode(filled)
        return GenomicArrayNode(ComputationNode(GenomicRunLengthArray.from_bedgraph,
                                                [interval_stream,
                                                 StreamNode(iter(genome_context.chrom_sizes.values()))]),
                                genome_context)

    def _get_intervals_from_data(self, name, data):
        if data.dtype == bool:
            return Interval([name] * len(data.starts),
                            data.starts, data.ends)[data.values]
        else:
            return BedGraph([name] * len(data.starts),
                            data.starts, data.ends, data.values)


class GenomicArrayGlobal(GenomicArray, np.lib.mixins.NDArrayOperatorsMixin):
    '''
    Class for holding the the enitre genomic array in memory
    '''

    def __init__(self, global_track: GenomicRunLengthArray, genome_context: GenomeContextBase):
        assert isinstance(global_track, GenomicRunLengthArray), global_track
        self._global_track = global_track
        self._genome_context = genome_context

    @property
    def dtype(self) -> np.dtype:
        return self._global_track.dtype

    def _index_boolean(self, idx):
        assert idx.dtype == bool
        assert isinstance(idx, GenomicArrayGlobal)
        return self._global_track[idx._global_track]

    def sum(self, axis=None) -> float:
        '''Sum the data in the array'''
        assert axis is None
        return self._global_track.sum(axis=None)

    def extract_chromsome(self, chromosome: Union[str, List[str]]) -> Union[
        GenomicRunLengthArray, RunLengthRaggedArray]:
        """Get the data on one or more chromosomes
        Parameters
        ----------
        chromosome : Union[str, List[str]]
            A chromosome name or a list of chromosome names

        Returns
        -------
        Union[GenomicRunLengthArray, RunLengthRaggedArray]
            A runlengtharray for a single chromosome or a RunLengthRaggedArray for a list of chromosomes
        """

        assert isinstance(chromosome, str)
        offset = self._genome_context.global_offset.get_offset([chromosome])[0]
        return self._global_track[offset:offset + self._genome_context.global_offset.get_size([chromosome])[0]]

    def __repr__(self) -> str:
        lines = []
        for name, _ in zip(self._genome_context.chrom_sizes, range(10)):
            lines.append(f'{name}: {self[name]}')
        if len(self._genome_context.chrom_sizes) > 10:
            lines.extend('...')
        return '\n'.join(lines)

    def to_dict(self) -> Dict[str, GenomicRunLengthArray]:
        """
        Convert the genomic array to a dict of arrays with chromosomes as keys

        Returns
        -------
        Dict[str, GenomicRunLengthArray]

        """
        go = self._genome_context.global_offset
        names = go.names()
        offsets = go.get_offset(names)
        sizes = go.get_size(names)
        return {name: self._global_track[offset:offset + size].to_array()
                for name, offset, size in zip(names, offsets, sizes)}

    def __array_ufunc__(self, ufunc: np.ufunc, method: str, *inputs, **kwargs):
        '''
        Handle numpy ufuncs on the genomic array

        Parameters
        ----------
        ufunc: np.ufunc
        method: str
            How to call the ufunc
        inputs: List
            Args for  the ufunc
        kwargs: Dict
            Additional arguments

        Returns
        -------
        GenomicArrayGlobal

        '''
        inputs = [(i._global_track if isinstance(i, GenomicArrayGlobal) else i) for i in inputs]
        r = self._global_track.__array_ufunc__(ufunc, method, *inputs, **kwargs)
        return self.__class__(r, self._genome_context)

    def __array_function__(self, func: Callable, types: List, args: List, kwargs: Dict):
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
        if func == np.sum:
            return self.sum(*args[1:], **kwargs)
        return NotImplemented

    def get_data(self) -> Union[Interval, BedGraph]:
        """The data of the array represented in BNPDataClass

        Returns
        -------
        Union[Interval, BedGraph]
            If boolean mask,return the intervals where True. Else a bedgraph containing the data
        """

        go = self._genome_context.global_offset
        names = go.names()
        starts = go.get_offset(names)
        stops = starts + go.get_size(names)
        intervals_list = []
        for name, start, stop in zip(names, starts, stops):
            assert isinstance(self._global_track, GenomicRunLengthArray)
            data = self._global_track[start:stop]
            assert isinstance(data, GenomicRunLengthArray)
            intervals_list.append(self._get_intervals_from_data(name, data))
        return np.concatenate(intervals_list)

    def extract_intervals(self, intervals: Union[Interval, 'GenomicIntervals'],
                          stranded: bool = False) -> RunLengthRaggedArray:
        """Extract the data contained in a set of intervals

        Parameters
        ----------
        intervals : Interval
        stranded : bool

        Returns
        -------
        RunLengthRaggedArray
        """
        global_intervals = self._genome_context.global_offset.from_local_interval(intervals)
        rle = self._global_track[global_intervals]
        if not stranded:
            return rle
        r = rle[:, ::-1]
        return np.where((intervals.strand.ravel() == '+')[:, np.newaxis],
                        rle, r)

    def extract_locations(self, locations: 'GenomicLocation', stranded=None) -> RunLengthRaggedArray:
        """Extract the data contained in a set of intervals

        Parameters
        ----------
        intervals : Interval
        stranded : bool

        Returns
        -------
        RunLengthRaggedArray
        """
        assert stranded is None, 'Stranded not implemented for locations'
        global_intervals = self._genome_context.global_offset.from_local_coordinates(locations.chromosome,
                                                                                     locations.position)
        return self._global_track[global_intervals]
        # if not stranded:
        #     return rle
        # r = rle[:, ::-1]
        # return np.where((intervals.strand.ravel() == '+')[:, np.newaxis],
        #                 rle, r)

    @classmethod
    def from_dict(cls, d: Dict[str, GenomicRunLengthArray], genome_context: GenomeContextBase) -> 'GenomicData':
        """Create genomic data from a dict of data with chromosomes as keys

        Parameters
        ----------
        d : Dict[str, GenomicRunLengthArray]

        Returns
        -------
        'GenomicData'
        """
        # chrom_sizes = {name: array.size for name, array in d.items()}
        array = np.concatenate(list(d.values()))
        return cls(array, genome_context)

    @classmethod
    def from_stream(cls, stream: Iterable[Tuple[str, GenomicRunLengthArray]],
                    genome_context: GenomeContextBase) -> 'GenomicData':
        '''
        Create a genomic array from a stream of data

        Parameters
        ----------
        stream: Iterable[Tuple[str, GenomicRunLengthArray]]
        genome_context: GenomeContextBase

        Returns
        -------
        GenomicData

        '''
        return cls.from_dict(dict(stream))


class GenomicArrayNode(GenomicArray, np.lib.mixins.NDArrayOperatorsMixin):
    '''
    Class for representing a genomic array, but only keeping data for one 
    chromosome in memory at a time
    '''

    def __str__(self):
        return 'GTN:' + str(self._run_length_node)

    __repr__ = __str__

    def __init__(self, run_length_node: ComputationNode, genome_context: GenomeContextBase):
        self._run_length_node = run_length_node
        self._chrom_name_node = StreamNode(iter(genome_context.chrom_sizes.keys()))
        self._genome_context = genome_context

    def __array_ufunc__(self, ufunc: np.ufunc, method: str, *inputs, **kwargs):
        args = [gtn._run_length_node if isinstance(gtn, GenomicArrayNode) else gtn for gtn in inputs]
        return self.__class__(ufunc(*args, **kwargs), self._genome_context)

    def get_data(self) -> ComputationNode:
        """Get the data from the array in a lazy BNPDataClass
        Returns
        -------
        ComputationNode
            A lazy compution node representing the data from the array
        """

        return ComputationNode(self._get_intervals_from_data, [self._chrom_name_node, self._run_length_node])

    def _index_boolean(self, idx):
        assert isinstance(idx, GenomicArrayNode)
        return ComputationNode(lambda a, i: a[i], [self._run_length_node, idx])

    def extract_intervals(self, intervals: Interval, stranded: bool = False) -> ComputationNode:
        """Extract the data on a set of intervals.

        Returns a `ComputationNode` representing the extracted data.

        Parameters
        ----------
        intervals : Interval
            The set of intervals to extract data from
        stranded : bool
            Whether or not to treat the intervals as stranded

        Returns
        -------
        ComputationNode:
            A computation node representing the raggedarray containing data from all the intervals
        """

        def stranded_func(ra: GenomicRunLengthArray, start: int, stop: int, strand: str):

            assert np.all(stop <= ra.size), (np.max(stop), ra.size)
            rle = ra[start:stop]
            r = rle[:, ::-1]
            return np.where((strand == '+')[:, np.newaxis],
                            rle, r)

        def unstranded_func(ra: GenomicRunLengthArray, start: int, stop: int):
            assert np.all(stop <= ra.size), (np.max(stop), ra.size)
            return ra[start:stop]

        assert self.genome_context.is_compatible(intervals.genome_context), (
        self.genome_context, intervals.genome_context)
        intervals = intervals.as_stream()
        if stranded:
            return ComputationNode(stranded_func, [self._run_length_node,
                                                   intervals.start,
                                                   intervals.stop,
                                                   intervals.strand])
        return ComputationNode(unstranded_func, [self._run_length_node, intervals.start, intervals.stop])

    def extract_chromsome(self, chromosome: Union[str, List[str]]) -> 'GenomicData':
        '''
        Extract the data for a chromosome
        Unimplemented
        '''
        assert False

    def from_stream(cls, stream: Iterable[Tuple[str, GenomicRunLengthArray]],
                    genome_context: GenomeContextBase) -> 'GenomicData':
        '''
        Create a genomic array from a stream of data

        Parameters
        ----------
        stream: Iterable[Tuple[str, GenomicRunLengthArray]]
        genome_context: GenomeContextBase

        Returns
        -------
        GenomicData

        '''
        stream_node = StreamNode((array for _, array in stream))
        return cls(stream_node, genome_context)

    @classmethod
    def from_dict(cls, d: Dict[str, GenomicRunLengthArray], genome_context: GenomeContextBase) -> 'GenomicData':
        """Create genomic data from a dict of data with chromosomes as keys

        Parameters
        ----------
        d : Dict[str, GenomicRunLengthArray]

        Returns
        -------
        'GenomicData'
        """
        # chrom_sizes = {name: array.size for name, array in d.items()}
        stream_node = StreamNode(iter(d.values()))
        return cls(stream_node, genome_context)

    def to_dict(self):
        '''Unimplemented'''
        assert False

    def sum(self, axis=None) -> float:
        '''Sum the data in the array'''
        return np.sum(self._run_length_node)

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
