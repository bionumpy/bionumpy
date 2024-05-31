from abc import ABC, abstractmethod, abstractproperty, abstractclassmethod
from collections import namedtuple

import numpy as np
from typing import List, Iterable, Tuple, Dict, Any, Optional

from .coordinate_mapping import map_locations, find_indices
from .. import EncodedArray
from ..bnpdataclass import BNPDataClass, replace, bnpdataclass
from .genomic_track import GenomicArray, GenomicArrayNode
from .genome_context_base import GenomeContextBase
from ..datatypes import Interval, Bed6, StrandedInterval, LocationEntry, StrandedLocationEntry
from ..arithmetics.intervals import get_pileup, merge_intervals, extend_to_size, clip, get_boolean_mask, RawInterval
from ..computation_graph import StreamNode, Node, ComputationNode, compute

from ..string_array import StringArray
from ..util.typing import EncodedArrayLike


class GenomicPlace:
    @property
    def genome_context(self):
        return self._genome_context

    @property
    @abstractmethod
    def get_location(self, where='start') -> 'GenomicLocation':
        return NotImplemented

    def get_data_field(self, field_name: str) -> Any:
        return NotImplemented

    def set_strand(self, strand: str):
        self._is_stranded = True
        self._strand = strand


class GenomicLocation(GenomicPlace):
    '''Class representing (possibliy stranded) locations in the genome'''

    @property
    @abstractmethod
    def chromosome(self):
        return NotImplemented

    @property
    @abstractmethod
    def position(self):
        return NotImplemented

    @property
    @abstractmethod
    def strand(self):
        return NotImplemented

    @property
    @abstractmethod
    def is_stranded(self):
        return NotImplemented

    @classmethod
    def from_fields(cls, genome_context: GenomeContextBase, chromosome: List[str], position: List[int],
                    strand: List[str] = None) -> 'GenomicLocation':
        """Create genomic location from a genome context and the needed fields (chromosome and position)

        Parameters
        ----------
        genome_context : GenomeContextBase
            Genome context object for the genome
        chromosome : List[str]
            List of chromosome
        position : List[int]:
            List of positions
        strand : List[str]
            Optional list of strand

        Returns
        -------
        'GenomicLocation'
        """
        is_stranded = strand is not None
        if is_stranded:
            data = StrandedLocationEntry(chromosome, position, strand)
        else:
            data = LocationEntry(chromosome, position)
        return GenomicLocationGlobal.from_data(data, genome_context, is_stranded=is_stranded)

    @classmethod
    def from_data(cls, data: BNPDataClass,
                  genome_context: GenomeContextBase,
                  is_stranded: bool = False,
                  chromosome_name: str = 'chromosome',
                  position_name: str = 'position',
                  strand_name: str = 'strand') -> 'GenomicLocation':
        """Create GenomicLocation object from a genome context and a bnpdataclass

        The field names for the chromosome, positions and strand can be specified

        Parameters
        ----------
        cls : 3
            4
        data : BNPDataClass
            The data containing the locations
        genome_context : GenomeContextBase
            Genome context for the genome
        is_stranded : bool
            Whether or not the locations should be stranded
        chromosome_name : str
            Name of the chromosome field in `data`
        position_name : str
            Name of the position field in the `data`
        strand_name : str
            Name if the `strand` field int the `data`

        Returns
        -------
        'GenomicLocation'
        """

        assert all(hasattr(data, name) for name in (chromosome_name, position_name))
        if is_stranded:
            assert hasattr(data, strand_name)
        return GenomicLocationGlobal(genome_context.mask_data(data), genome_context, is_stranded,
                                     {'chromosome': chromosome_name,
                                      'position': position_name,
                                      'strand': strand_name})


class GenomicLocationGlobal(GenomicLocation, ABC):
    ''' Class for genomic locations that are kept entirely in memory'''

    def __init__(self, locations: BNPDataClass, genome_context: GenomeContextBase, is_stranded: bool,
                 field_dict: Dict[str, str]):
        self._locations = locations
        self._genome_context = genome_context
        self._is_stranded = is_stranded
        self._field_dict = field_dict

    def __repr__(self):
        return f'Genomic Locations on {self._genome_context}:\n{self._locations.astype(LocationEntry)}'

    @property
    def data(self) -> BNPDataClass:
        '''
        Return the underlying data as a bnpdataclass

        Returns
        -------
        BNPDataClass
        '''
        return self._locations

    def __replace__(self, **kwargs):
        '''
        Replace fields in the locations, used internally by bnp.replace
        '''

        kwargs = {self._field_dict[kw]: value for kw, value in kwargs.items()}
        return self.__class__(replace(self._locations, **kwargs), self._genome_context, self._is_stranded,
                              self._field_dict)

    @property
    def chromosome(self) -> StringArray:
        '''
        Return the chromosome of the locations

        Returns
        -------
        StringArray

        '''
        return getattr(self._locations, self._field_dict['chromosome'])

    @property
    def position(self) -> np.ndarray:
        """
        Return the position of the locations

        Returns
        -------
        np.ndarray

        """
        return getattr(self._locations, self._field_dict['position'])

    @property
    def strand(self) -> EncodedArrayLike:
        """
        The strand of the locations

        Returns
        -------
        EncodedArrayLike

        """
        if not self.is_stranded():
            raise ValueError('Unstranded position has not strand')
        return getattr(self._locations, self._field_dict['strand'])

    def is_stranded(self) -> bool:
        '''
        Return whether the locations are stranded

        Returns
        -------
        bool

        '''
        return self._is_stranded

    def get_windows(self, flank: Optional[int] = None, window_size: Optional[int] = None) -> 'GenomicIntervals':
        """Create windows around the locations. 

        `Flank specifies the flank on either side of the location. The full windows
        will thus be `flank*2+1` wide

        Parameters
        ----------
        flank : int
            Flank on either side of the location

        Returns
        -------
        GenomicIntervals
            Window intervals

        """
        if flank is not None:
            assert window_size is None
            l_flank = flank
            r_flank = flank + 1
        else:
            assert window_size is not None
            l_flank = window_size // 2
            r_flank = window_size // 2 + window_size % 2
        if self.is_stranded():
            intervals = StrandedInterval(self.chromosome, self.position - l_flank,
                                         self.position + r_flank, self.strand)
        else:
            intervals = Interval(self.chromosome, self.position - l_flank,
                                 self.position + r_flank)
        return GenomicIntervalsFull(intervals, self._genome_context,
                                    is_stranded=self.is_stranded()).clip()

    def sorted(self) -> GenomicLocation:
        """Return a sorted version of the locations

        Sorted according the chromosome order in the `GenomeContext`

        Returns
        -------
        GenomicLocation
            Sorted locations

        """
        return self[np.lexsort([self.position, self.chromosome.raw()])]

    def __getitem__(self, idx: Any)-> 'GenomicLocation':
        '''
        Get a subset of the locations

        Parameters
        ----------
        idx: numpy index like

        Returns
        -------
        GenomicLocation

        '''
        return self.__class__(self._locations[idx], self._genome_context, self._is_stranded, self._field_dict)

    def get_data_field(self, field_name: str) -> Any:
        '''
        Get a field from the undelying data

        Parameters
        ----------
        field_name: str

        Returns
        -------
        Any
            The field from the underlying data

        '''
        return getattr(self._locations, field_name)


class GenomicLocationStreamed(GenomicLocation, Node):
    '''
    Class for representing intervals that are grouped by chromosome, and where only intervals
    for one chromosome at the time is kept in memory
    '''

    is_stream = True

    def _get_chrom_size(self, intervals: Interval):
        return self._genome_context.chrom_sizes[intervals.chromosome]

    def __str__(self):
        return 'GLS:' + str(self._data_node)

    def __repr__(self):
        return 'GLS:' + str(self._data_node)

    def __init__(self, data_node: Node, genome_context: GenomeContextBase, is_stranded=False,
                 field_dict: Dict[str, str] = None):
        if field_dict is None:
            field_dict = {name: name for name in ['chromosome', 'positions', 'strand']}
        self._genome_context = genome_context
        self._chromosome = ComputationNode(getattr, [data_node, field_dict['chromosome']])
        self._position = ComputationNode(getattr, [data_node, field_dict['position']])
        if is_stranded:
            self._strand = ComputationNode(getattr, [data_node, field_dict['strand']])

        self._chrom_size_node = StreamNode(iter(self._genome_context.chrom_sizes.values()))
        self._data_node = data_node
        self._is_stranded = is_stranded

    def is_stranded(self) -> bool:
        return self._is_stranded

    def sorted(self) -> 'GenomicLocationStreamed':
        return NotImplemented

    @property
    def position(self) -> ComputationNode:
        '''
        Return a computation node for the position of the locations.
        Can be realized by calling bnp.compute on them, or can be aggregated across the stream with
        numpy aggregations
        -------
        ComputationNode
        '''
        return self._position

    @property
    def chromosome(self) -> ComputationNode:
        '''
        Return a computation node for the chromosome of the locations.
        Can be realized by calling bnp.compute on them, or can be aggregated across the stream with
        numpy aggregations

        Returns
        -------
        ComputationNode

        '''
        return self._chromosome

    def get_data_field(self, field_name: str) -> ComputationNode:
        '''
        Get a field from the undelying data returned as a computation node

        Parameters
        ----------
        field_name: str

        Returns
        -------
        ComputationNode

        '''
        return ComputationNode(getattr, [self._data_node, field_name])

    @property
    def strand(self) -> ComputationNode:
        '''
        Return a computation node for the strand of the locations.
        Can be realized by calling bnp.compute on them, or can be aggregated across the stream with
        numpy aggregations

        Returns
        -------
        ComputationNode

        '''

        if not self.is_stranded():
            raise ValueError('Strand not supported on unstranded intervals')
        return self._strand

    def __getitem__(self, item):
        '''
        Get a subset of the locations, returned as a CompuationNode

        Parameters
        ----------
        item

        Returns
        -------

        '''
        return self.__class__(ComputationNode(lambda x, i: x[i], [self._data_node, item]), self._genome_context)

    def get_windows(self, flank: int = None, window_size: Optional[int] = None) -> 'GenomicIntervals':
        """Create windows around the locations. 

        `Flank specifies the flank on either side of the location. The full windows
        will thus be `flank*2+1` wide

        Parameters
        ----------
        window_size
        flank : int
            Flank on either side of the location

        Returns
        -------
        GenomicIntervals
            Window intervals

        """
        if flank is not None:
            assert window_size is None
            l_flank = flank
            r_flank = flank + 1
        else:
            assert window_size is not None
            l_flank = window_size // 2
            r_flank = window_size // 2 + window_size % 2
        if self.is_stranded():
            intervals = ComputationNode(StrandedInterval,
                                        [self.chromosome,
                                         self.position - l_flank,
                                         self.position + r_flank, self.strand])
        else:
            intervals = ComputationNode(
                Interval, [self.chromosome, self.position - l_flank,
                           self.position + r_flank])
        return GenomicIntervalsStreamed(intervals, self._genome_context,
                                        is_stranded=self.is_stranded()).clip()

    def _get_buffer(self, i):
        return GenomicLocationGlobal(
            LocationEntry(self.chromosome._get_buffer(i),
                          self.position._get_buffer(i)),
            self._genome_context)


class GenomicIntervals(GenomicPlace):
    ''' Class for representing intervals on a genome'''

    @property
    @abstractmethod
    def start(self):
        return NotImplemented

    @property
    @abstractmethod
    def stop(self):
        return NotImplemented

    @property
    @abstractmethod
    def chromosome(self):
        return NotImplemented

    @property
    @abstractmethod
    def strand(self):
        return NotImplemented

    @abstractmethod
    def is_stranded(self):
        return NotImplemented

    @abstractmethod
    def get_location(self, where: str = 'start') -> GenomicLocation:
        return NotImplemented

    @abstractmethod
    def extended_to_size(self, size: int) -> 'GenomicIntervals':
        """Extend intervals along strand to reach the given size

        Parameters
        ----------
        size : int

        Returns
        -------
        'GenomicIntervals'
        """
        return NotImplemented

    @abstractmethod
    def merged(self, distance: int = 0) -> 'GenomicIntervals':
        """Merge intervals that overlap or lie within distance of eachother

        Parameters
        ----------
        distance : int

        Returns
        -------
        'GenomicIntervals'
            4
        """

        return NotImplemented

    @abstractmethod
    def get_mask(self) -> GenomicArray:
        """Return a boolean mask of areas covered by any interval

        Returns
        -------
        GenomicArray
            Genomic mask
        """

        return NotImplemented

    @abstractmethod
    def get_pileup(self) -> GenomicArray:
        """Return a genmic track of counting the number of intervals covering each bp

        Returns
        -------
        GenomicArray
            Pileup track
        """
        return NotImplemented

    @classmethod
    def from_track(cls, track: GenomicArray) -> 'GenomicIntervals':
        """Return intervals of contigous areas of nonzero values of track

        Parameters
        ----------
        track : GenomicArray

        Returns
        -------
        'GenomicIntervals'
        """
        if isinstance(track, GenomicArrayNode):
            return GenomicIntervalsStreamed(track.get_data(), track.genome_context)
        return GenomicIntervalsFull(track.get_data(), track.genome_context)

    @classmethod
    def from_fields(cls, genome_context: GenomeContextBase, chromosome: StringArray, start: np.ndarray, stop: np.ndarray, strand: Optional[EncodedArray]=None) -> 'GenomicIntervals':
        '''
        Create genomic intervals from fields

        Parameters
        ----------
        genome_context: GenomeContextBase
        chromosome: StringArray
        start: np.ndarray
        stop: np.ndarray
        strand: EncodedArray

        Returns
        -------
        GenomicIntervals
        '''

        is_stranded = strand is not None
        if is_stranded:
            intervals = Bed6(chromosome, start, stop, ['.'] * len(start),
                             np.zeros_like(start), strand)
        else:
            intervals = Interval(chromosome, start, stop)
        return cls.from_intervals(intervals, genome_context, is_stranded=is_stranded)

    @classmethod
    def from_intervals(cls, intervals: Interval, genome_context: GenomeContextBase, is_stranded: Optional[bool]=False) -> 'GenomicIntervals':
        """Create genomic intervals from interval entries and genome info

        Parameters
        ----------
        intervals : Interval
        genome_context : GenomeContextBase
        is_stranded : bool

        Returns
        -------
        'GenomicIntervals'

        """
        if isinstance(intervals, Interval):  # TODO check is node
            return GenomicIntervalsFull(genome_context.mask_data(intervals), genome_context, is_stranded)
        else:
            return cls.from_interval_stream(intervals, genome_context, is_stranded)

    @classmethod
    def from_interval_stream(cls, interval_stream: Iterable[Interval], genome_context: GenomeContextBase,
                             is_stranded=False) -> 'GenomicIntervals':
        """Create streamed genomic intervals from a stream of intervals and genome info

        Parameters
        ----------
        interval_stream : Iterable[Interval]
        genome_context : GenomeContextBase
        is_stranded : bool

        Returns
        -------
        'GenomicIntervals'
        """

        interval_stream = genome_context.iter_chromosomes(
            interval_stream, StrandedInterval if is_stranded else Interval)
        return GenomicIntervalsStreamed(StreamNode(interval_stream), genome_context, is_stranded=is_stranded)

    @abstractmethod
    def clip(self) -> 'GenomicIntervals':
        """Clip the intervals so that they are contained in the genome

        Returns
        -------
        'GenomicIntervals'
            Clipped intervals
        """
        return NotImplemented

    def compute(self):
        return NotImplemented


class GenomicIntervalsFull(GenomicIntervals):
    ''' Class for holding a set of intervals in memory'''

    is_stream = False

    def __init__(self, intervals: Interval, genome_context: GenomeContextBase, is_stranded=False):
        self._intervals = intervals
        self._is_stranded = is_stranded
        self._genome_context = genome_context

    @property
    def data(self) -> BNPDataClass:
        '''Return the underlying data as a bnpdataclass'''
        return self._intervals

    def __array_function__(self, func: callable, types: List, args: List, kwargs: Dict):
        if func == np.concatenate:
            return self.__class__(np.concatenate([obj._intervals for obj in args[0]]), self._genome_context,
                                  self._is_stranded)

        return NotImplemented

    def __repr__(self):
        return f'Genomic Intervals on {self._genome_context}:\n{self._intervals.astype(Interval)}'

    def get_data(self) -> BNPDataClass:
        """Return the underlying data for the intervals

        Returns
        -------
        BNPDataClass
            The data for the intervals

        """
        return self._intervals

    def __len__(self) -> int:
        return len(self._intervals)

    def map_locations(self, locations: LocationEntry):
        '''
        Map locations to the intervals. The locations should be in the same coordinate system as the intervals
        The new locations will be in the coordinate system of the intervals

        Parameters
        ----------
        locations: LocationEntry

        Returns
        -------
        LocationEntry

        '''
        go = self._genome_context.global_offset.from_local_interval(self._intervals)
        global_positions = self._genome_context.global_offset.from_local_coordinates(locations.chromosome,
                                                                                     locations.position)

        location_indices, interval_indices = find_indices(global_positions, go)
        new_entries = locations[location_indices]
        names = self._intervals.name if hasattr(self._intervals, 'name') else StringArray(
            np.arange(len(self._intervals)).astype('S'))
        return replace(new_entries, chromosome=names[interval_indices],
                       position=new_entries.position - self.start[interval_indices])

        return map_locations(replace(locations, position=global_positions), go)

    def sorted(self) -> 'GenomicIntervals':
        """Return the intervals sorted according to `genome_context`

        Returns
        -------
        'GenomicIntervals'

        """
        args = np.lexsort([self.stop, self.start, self.chromosome.raw()])
        return self[args]

    def __getitem__(self, idx):
        return self.__class__(self._intervals[idx], self._genome_context, self._is_stranded)

    def get_location(self, where: str = 'start') -> GenomicLocation:
        """Get the genomic location of eitert 'start', 'stop' or 'center' of the intervals

        Parameters
        ----------
        where : str
            'start', 'stop' or 'center'

        Returns
        -------
        GenomicLocation
        """

        if where in ('start', 'stop'):
            if not self.is_stranded():
                data = self._intervals
            else:
                location = np.where(self.strand == ('+' if where == 'start' else '-'),
                                    self.start,
                                    self.stop - 1)
                data = replace(self._intervals, start=location)
        else:
            assert where == 'center'
            location = (self.start + self.stop) // 2
            data = replace(self._intervals, start=location)
        return GenomicLocationGlobal.from_data(
            data, self._genome_context, is_stranded=self.is_stranded(),
            position_name='start')

    @property
    def start(self) -> int:
        return self._intervals.start

    @property
    def stop(self) -> int:
        return self._intervals.stop

    @property
    def strand(self) -> str:
        if not self.is_stranded():
            raise ValueError('Unstranded interval has not strand')
        return self._intervals.strand

    def get_data_field(self, field_name: str):
        return getattr(self._intervals, field_name)

    @property
    def chromosome(self) -> str:
        return self._intervals.chromosome

    def extended_to_size(self, size: int) -> GenomicIntervals:
        """Extend intervals along strand to reach the given size

        Parameters
        ----------
        size : int

        Returns
        -------
        'GenomicIntervals'
        """
        chrom_sizes = self._genome_context.global_offset.get_size(self._intervals.chromosome)
        return self.from_intervals(extend_to_size(self._intervals, size, chrom_sizes),
                                   self._genome_context)

    def merged(self, distance: int = 0) -> GenomicIntervals:
        """Merge intervals that overlap or lie within distance of eachother

        Parameters
        ----------
        distance : int

        Returns
        -------
        'GenomicIntervals'
        """

        if distance > 0:
            stream = self.as_stream()
            return stream.merged(distance).compute()

        assert distance == 0, 'Distance might cross chromosome boundries so is not supported with current implementation'
        go = self._genome_context.global_offset
        global_intervals = go.from_local_interval(self._intervals)
        global_merged = merge_intervals(global_intervals, distance)
        return self.from_intervals(
            self._global_offset.to_local_interval(global_merged), self._genome_context)

    def get_pileup(self) -> GenomicArray:
        """Return a genmic array of counting the number of intervals covering each bp

        Returns
        -------
        GenomicArray
            Pileup track
        """
        go = self._genome_context.global_offset.from_local_interval(self._intervals)
        return GenomicArray.from_global_data(
            get_pileup(go, self._genome_context.size),
            self._genome_context)

    def get_mask(self) -> GenomicArray:
        """Return a boolean mask of areas covered by any interval

        Returns
        -------
        GenomicArray
            Genomic mask
        """
        I = RawInterval
        starts, stops = self._genome_context.global_offset.start_ends_from_intervals(self._intervals)
        global_mask = get_boolean_mask(I(starts, stops), self._genome_context.size)
        return GenomicArray.from_global_data(global_mask, self._genome_context)

    def clip(self) -> 'GenomicIntervalsFull':
        """Clip the intervals so that they are contained in the genome

        Returns
        -------
        'GenomicIntervals'
            Clipped intervals
        """
        chrom_sizes = self._genome_context.global_offset.get_size(self._intervals.chromosome)
        return replace(self,
                       start=np.maximum(0, self.start),
                       stop=np.minimum(chrom_sizes, self.stop))

    def __replace__(self, **kwargs):
        return self.__class__(replace(self._intervals, **kwargs), self._genome_context, self._is_stranded)

    def compute(self):
        return self

    def as_stream(self):
        interval_class = StrandedInterval if self._is_stranded else Interval
        filled = self.genome_context.iter_chromosomes(self._intervals, interval_class)
        return GenomicIntervalsStreamed(
            StreamNode(filled),
            self._genome_context, self._is_stranded)

    def get_sorted_stream(self):
        sorted_intervals = self.sorted()
        return self.from_interval_stream(iter([sorted_intervals]))

    def is_stranded(self):
        return self._is_stranded


class GenomicIntervalsStreamed(GenomicIntervals, Node):
    '''
    Class for representing intervals that are grouped by chromosome, and where only intervals
    for one chromosome at the time is kept in memory
    '''

    is_stream = True

    def _get_chrom_size(self, intervals: Interval):
        return self._genome_context.chrom_sizes[intervals.chromosome]

    def __str__(self):
        return 'GIS:' + str(self._intervals_node)

    def __repr__(self):
        return 'GIS:' + str(self._intervals_node)

    def __init__(self, intervals_node: Node, genome_context: GenomeContextBase, is_stranded=False):
        self._genome_context = genome_context
        self._start = ComputationNode(getattr, [intervals_node, 'start'])
        self._stop = ComputationNode(getattr, [intervals_node, 'stop'])
        if is_stranded:
            self._strand = ComputationNode(getattr, [intervals_node, 'strand'])
        self._chromosome = ComputationNode(getattr, [intervals_node, 'chromosome'])
        self._chrom_size_node = StreamNode(iter(self._genome_context.chrom_sizes.values()))
        self._intervals_node = intervals_node
        self._is_stranded = is_stranded

    def is_stranded(self):
        return self._is_stranded

    def sorted(self):
        return NotImplemented

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def chromosome(self):
        return self._chromosome

    def get_data_field(self, field_name: str):
        return ComputationNode(getattr, [self._intervals_node, 'chromosome'])

    @property
    def strand(self):
        if not self.is_stranded():
            raise ValueError('Strand not supported on unstranded intervals')
        return self._strand

    def __getitem__(self, item):
        return self.__class__(ComputationNode(lambda x, i: x[i], [self._intervals_node, item]), self._genome_context)

    def extended_to_size(self, size: int) -> GenomicIntervals:
        """Extend intervals along strand to reach the given size

        Parameters
        ----------
        size : int

        Returns
        -------
        'GenomicIntervals'
        """
        return self.__class__(
            ComputationNode(extend_to_size, [self._intervals_node, size, self._chrom_size_node]),
            self._genome_context)

    def merged(self, distance: int = 0) -> GenomicIntervals:
        """Merge intervals that overlap or lie within distance of eachother

        Parameters
        ----------
        distance : int

        Returns
        -------
        'GenomicIntervals'
            4
        """

        return self.__class__(ComputationNode(merge_intervals, [self._intervals_node, distance]), self._genome_context)

    def get_pileup(self) -> GenomicArray:
        """Create a GenomicTrack of how many intervals covers each position in the genome

        Parameters
        ----------
        intervals : Interval

        Returns
        -------
        GenomicArray
        """
        return GenomicArrayNode(ComputationNode(get_pileup, [self._intervals_node, self._chrom_size_node]),
                                self._genome_context)

    def get_mask(self) -> GenomicArray:
        return GenomicArrayNode(ComputationNode(get_boolean_mask, [self._intervals_node, self._chrom_size_node]),
                                self._genome_context)

    def clip(self) -> 'GenomicIntervals':
        return self.__class__(ComputationNode(clip, [self._intervals_node, self._chrom_size_node]),
                              self._genome_context)

    def __replace__(self, **kwargs):
        return self.__class__(ComputationNode(replace, [self._intervals_node], kwargs), self._genome_context)
        # return self.__class__(dataclasses.replace(self._intervals, **kwargs), self._genome_context)

    def compute(self):
        chromosome, start, stop = compute((self.chromosome, self.start, self.stop))
        return GenomicIntervalsFull(Interval(chromosome, start, stop), self._genome_context)

    def _get_buffer(self, i):
        return GenomicIntervalsFull(Interval(self.chromosome._get_buffer(i),
                                             self.start._get_buffer(i),
                                             self.stop._get_buffer(i)),
                                    self._genome_context)

    def as_stream(self):
        return self

    def get_location(self, where: str = 'start') -> GenomicLocation:
        """Get the genomic location of eitert 'start', 'stop' or 'center' of the intervals

        Parameters
        ----------
        where : str
            'start', 'stop' or 'center'

        Returns
        -------
        GenomicLocation
        """
        assert where == 'start' and not self.is_stranded()
        return GenomicLocationStreamed(self._intervals_node,
                                       self._genome_context,
                                       False,
                                       {'chromosome': 'chromosome',
                                        'position': 'start',
                                        'strand': 'strand'})
