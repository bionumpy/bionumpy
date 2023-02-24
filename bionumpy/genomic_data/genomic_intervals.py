from abc import ABC, abstractmethod, abstractproperty, abstractclassmethod
import numpy as np
from typing import List, Iterable, Tuple, Dict
from ..bnpdataclass import BNPDataClass, replace, bnpdataclass
from ..encodings import StrandEncoding
from .genomic_track import GenomicArray, GenomicArrayNode
from .genome_context_base import GenomeContextBase
from ..datatypes import Interval, Bed6, StrandedInterval
from ..arithmetics.intervals import get_pileup, merge_intervals, extend_to_size, clip, get_boolean_mask
from ..computation_graph import StreamNode, Node, ComputationNode, compute
import dataclasses


@bnpdataclass
class LocationEntry:
    chromosome: str
    position: int


@bnpdataclass
class StrandedLocationEntry(LocationEntry):
    strand: StrandEncoding


class GenomicPlace:
    @property
    def genome_context(self):
        return self._genome_context

    @abstractproperty
    def get_location(self, where='start'):
        return NotImplemented


class GenomicLocation(GenomicPlace):
    @abstractproperty
    def chromosome(self):
        return NotImplemented

    @abstractproperty
    def position(self):
        return NotImplemented

    @abstractproperty
    def strand(self):
        return NotImplemented

    @abstractmethod
    def is_stranded(self):
        return NotImplemented

    @classmethod
    def from_fields(cls, genome_context, chromosome, position, strand=None):
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
                  strand_name: str = 'strand'):
        assert all(hasattr(data, name) for name in (chromosome_name, position_name))
        if is_stranded:
            assert hasattr(data, strand_name)
        return GenomicLocationGlobal(data, genome_context, is_stranded,
                                     {'chromosome': chromosome_name,
                                      'position': position_name,
                                      'strand': strand_name})


class GenomicLocationGlobal(GenomicLocation):
    def __init__(self, locations, genome_context: GenomeContextBase, is_stranded, field_dict):
        self._locations = locations
        self._genome_context = genome_context
        self._is_stranded = is_stranded
        self._field_dict = field_dict

    @property
    def chromosome(self):
        return getattr(self._locations, self._field_dict['chromosome'])

    @property
    def position(self):
        return getattr(self._locations, self._field_dict['position'])

    @property
    def strand(self):
        if not self.is_stranded():
            raise ValueError('Unstranded position has not strand')
        return getattr(self._locations, self._field_dict['strand'])

    def is_stranded(self):
        return self._is_stranded

    def get_windows(self, flank):
        if self.is_stranded():
            intervals = StrandedInterval(self.chromosome, self.position-flank,
                                         self.position+flank, self.strand)
        else:
            intervals = Interval(self.chromosome, self.position-flank,
                                 self.position+flank)
        return GenomicIntervalsFull(intervals, self._genome_context, is_stranded=self.is_stranded()).clip()
    # self._genome_context.geometry.clip(intervals), self._genome_context,
    # is_stranded=self.is_stranded())


class GenomicIntervals(GenomicPlace):
    @abstractproperty
    def start(self):
        return NotImplemented

    @abstractproperty
    def stop(self):
        return NotImplemented

    @abstractproperty
    def chromosome(self):
        return NotImplemented

    @abstractproperty
    def strand(self):
        return NotImplemented

    @abstractmethod
    def is_stranded(self):
        return NotImplemented

    @abstractmethod
    def get_location(self, where='start'):
        return NotImplemented

    @abstractmethod
    def extended_to_size(self, size: int) -> 'GenomicIntervals':
        """Extend intervals along strand to rach the given size

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
            return GenomicIntervalsStreamed(track.get_data(), track._genome_context)
        return GenomicIntervalsFull(track.get_data(), track._genome_context)

    @classmethod
    def from_fields(cls, genome_context: GenomeContextBase, chromosome, start, stop, strand=None):
        is_stranded = strand is not None
        if is_stranded:
            intervals = Bed6(chromosome, start, stop, ['.']*len(start),
                             np.zeros_like(start), strand)
        else:
            intervals = Interval(chromosome, start, stop)
        return cls.from_intervals(intervals, genome_context, is_stranded=is_stranded)

    @classmethod
    def from_intervals(cls, intervals: Interval, genome_context: GenomeContextBase, is_stranded=False):
        """Create genomic intervals from interval entries and genome info

        Parameters
        ----------
        intervals : Interval
        chrom_sizes : Dict[str, int]
        """
        if isinstance(intervals, Interval): #TODO check is node
            return GenomicIntervalsFull(genome_context.mask_data(intervals), genome_context, is_stranded)
        else:
            return cls.from_interval_stream(intervals, genome_context, is_stranded)

    @classmethod
    def from_interval_stream(cls, interval_stream: Iterable[Interval], genome_context: GenomeContextBase, is_stranded=False):
        """Create streamed genomic intervals from a stream of intervals and genome info

        Parameters
        ----------
        interval_stream : Iterable[Interval]
        chrom_sizes : Dict[str, int]
        """
        
        interval_stream = genome_context.iter_chromosomes(
            interval_stream, StrandedInterval if is_stranded else Interval)
        return GenomicIntervalsStreamed(StreamNode(interval_stream), genome_context, is_stranded=is_stranded)

    @abstractmethod
    def clip(self) -> 'GenomicIntervals':
        return NotImplemented

    def compute(self):
        return NotImplemented


class GenomicIntervalsFull(GenomicIntervals):
    is_stream=False

    def __init__(self, intervals: Interval, genome_context: GenomeContextBase, is_stranded=False):
        self._intervals = intervals
        self._is_stranded = is_stranded
        self._genome_context = genome_context

    def __repr__(self):
        return f'Genomic Intervals on {self._genome_context}:\n{self._intervals.astype(Interval)}'

    def get_data(self):
        return self._intervals

    def __len__(self):
        return len(self._intervals)

    def sorted(self):
        args = np.lexsort([self.stop, self.start, self.chromosome.raw()])
        return self[args]

    def __getitem__(self, idx):
        return self.__class__(self._intervals[idx], self._genome_context, self._is_stranded)

    def get_location(self, where='start'):
        if where in ('start', 'stop'):
            if not self.is_stranded():
                data = self._intervals
            else:
                location = np.where(self.strand==('+' if where=='start' else '-'),
                                    self.start,
                                    self.stop-1)
                data = replace(self._intervals, start=location)
        else:
            assert where == 'center'
            location = (self.start+self.stop)//2
            data = replace(self._intervals, start=location)
        return GenomicLocationGlobal.from_data(
            data, self._genome_context, is_stranded=self.is_stranded(),
            position_name='start')


    @property
    def start(self):
        return self._intervals.start

    @property
    def stop(self):
        return self._intervals.stop

    @property
    def strand(self):
        if not self.is_stranded():
            raise ValueError('Unstranded interval has not strand')
        return self._intervals.strand

    @property
    def chromosome(self):
        return self._intervals.chromosome

    def extended_to_size(self, size: int) -> GenomicIntervals:
        chrom_sizes = self._genome_context.global_offset.get_size(self._intervals.chromosome)
        return self.from_intervals(extend_to_size(self._intervals, size, chrom_sizes), 
                                   self._genome_context)

    def merged(self, distance: int = 0) -> GenomicIntervals:
        assert distance == 0, 'Distance might cross chromosome boundries so is not supported with current implementation'
        go = self._genome_context.global_offset
        global_intervals = go.from_local_interval(self._intervals)
        global_merged = merge_intervals(global_intervals, distance)
        return self.from_intervals(
            self._global_offset.to_local_interval(global_merged), self._genome_context)

    def get_pileup(self) -> GenomicArray:
        go = self._genome_context.global_offset.from_local_interval(self._intervals)
        return GenomicArray.from_global_data(
            get_pileup(go, self._genome_context.size),
            self._genome_context)

    def get_mask(self) -> GenomicArray:
        go = self._genome_context.global_offset.from_local_interval(self._intervals)
        global_mask = get_boolean_mask(go, self._genome_context.size)
        return GenomicArray.from_global_data(global_mask, self._genome_context)

    def clip(self) -> 'GenomicIntervalsFull':
        chrom_sizes = self._genome_context.global_offset.get_size(self._intervals.chromosome)
        return replace(self,
                       start=np.maximum(0, self.start),
                       stop=np.minimum(chrom_sizes, self.stop))

    def __replace__(self, **kwargs):
        return self.__class__(dataclasses.replace(self._intervals, **kwargs), self._genome_context)

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


class GenomicIntervalsStreamed(GenomicIntervals):
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

    @property
    def strand(self):
        if not self.is_stranded():
            raise ValueError('Strand not supported on unstranded intervals')
        return self._strand

    def __getitem__(self, item):
        return self.__class__(ComputationNode(lambda x, i: x[i], [self._intervals_node, item]), self._genome_context)

    def extended_to_size(self, size: int) -> GenomicIntervals:
        return self.__class__(
            ComputationNode(extend_to_size, [self._intervals_node, size, self._chrom_size_node]),
            self._genome_context)

    def merged(self, distance: int = 0) -> GenomicIntervals:
        return self.__class__(ComputationNode(merge_intervals, [self._intervals_node]), self._genome_context)

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

    def clip(self) -> 'GenomicIntervalsFull':
        return self.__class__(ComputationNode(clip, [self._intervals_node, self._chrom_size_node]), self._genome_context)


    def __replace__(self, **kwargs):
        return self.__class__(ComputationNode(dataclasses.replace, [self._intervals_node], kwargs), self._genome_context)
        return self.__class__(dataclasses.replace(self._intervals, **kwargs), self._genome_context)

    def compute(self):
        chromosome, start, stop = compute(self.chromosome, self.start, self.stop)
        return GenomicIntervalsFull(Interval(chromosome, start, stop), self._genome_context)

    def as_stream(self):
        return self
