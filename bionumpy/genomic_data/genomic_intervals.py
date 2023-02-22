from abc import ABC, abstractmethod, abstractproperty, abstractclassmethod
import numpy as np
from typing import List, Iterable, Tuple, Dict
from ..bnpdataclass import BNPDataClass, replace, bnpdataclass
from ..encodings import StrandEncoding
from .genomic_track import GenomicArray, GenomicArrayNode, fill_grouped, GenomeContext
from ..datatypes import Interval, Bed6, StrandedInterval
from ..arithmetics.intervals import get_pileup, merge_intervals, extend_to_size, clip, get_boolean_mask
from ..streams import groupby
from ..computation_graph import StreamNode, Node, ComputationNode, compute
from .geometry import Geometry
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
        return GenomeContext.from_dict(self._chrom_sizes)

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
    def from_fields(cls, chrom_sizes, chromosome, position, strand=None):
        is_stranded = strand is not None
        if is_stranded:
            data = StrandedLocationEntry(chromosome, position, strand)
        else:
            data = LocationEntry(chromosome, position)
        return GenomicLocationGlobal.from_data(data, chrom_sizes, is_stranded=is_stranded)

    @classmethod
    def from_data(cls, data: BNPDataClass,
                  chrom_sizes: Dict[str, int],
                  is_stranded: bool = False,
                  chromosome_name: str = 'chromosome',
                  position_name: str = 'position',
                  strand_name: str = 'strand'):
        assert all(hasattr(data, name) for name in (chromosome_name, position_name))
        if is_stranded:
            assert hasattr(data, strand_name)
        return GenomicLocationGlobal(data, chrom_sizes, is_stranded,
                                     {'chromosome': chromosome_name,
                                      'position': position_name,
                                      'strand': strand_name})


class GenomicLocationGlobal(GenomicLocation):
    def __init__(self, locations, chrom_sizes, is_stranded, field_dict):
        self._locations = locations
        self._genome_context = GenomeContext.from_dict(chrom_sizes)
        self._chrom_sizes = chrom_sizes
        self._geometry = Geometry(self._genome_context.chrom_sizes)
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
        return GenomicIntervalsFull(
            self._geometry.clip(intervals), self._chrom_sizes,
            is_stranded=self.is_stranded())


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
    def from_fields(cls, chrom_sizes, chromosome, start, stop, strand=None):
        is_stranded = strand is not None
        if is_stranded:
            intervals = Bed6(chromosome, start, stop, ['.']*len(start),
                             np.zeros_like(start), strand)
        else:
            intervals = Interval(chromosome, start, stop)
        return cls.from_intervals(intervals, chrom_sizes, is_stranded=is_stranded)

    @classmethod
    def from_intervals(cls, intervals: Interval, chrom_sizes: Dict[str, int], is_stranded=False):
        """Create genomic intervals from interval entries and genome info

        Parameters
        ----------
        intervals : Interval
        chrom_sizes : Dict[str, int]
        """
        if isinstance(intervals, Interval): #TODO check is node
            return GenomicIntervalsFull(intervals, chrom_sizes, is_stranded)
        else:
            return cls.from_interval_stream(intervals, chrom_sizes, is_stranded)

    @classmethod
    def from_interval_stream(cls, interval_stream: Iterable[Interval], chrom_sizes: Dict[str, int], is_stranded=False):
        """Create streamed genomic intervals from a stream of intervals and genome info

        Parameters
        ----------
        interval_stream : Iterable[Interval]
        chrom_sizes : Dict[str, int]
        """
        
        # filled = fill_grouped(groupby(bedgraph, 'chromosome'), chrom_sizes.keys(), BedGraph)
        interval_stream = GenomeContext.from_dict(chrom_sizes).iter_chromosomes(interval_stream, StrandedInterval if is_stranded else Interval)
        # grouped = groupby(interval_stream, 'chromosome')
        # interval_stream = StreamNode(fill_grouped(grouped, chrom_sizes.keys(), StrandedInterval if is_stranded else Interval))
        # pair[1] for pair in interval_stream)
        return GenomicIntervalsStreamed(StreamNode(interval_stream), chrom_sizes, is_stranded=is_stranded)

    @abstractmethod
    def clip(self) -> 'GenomicIntervals':
        return NotImplemented

    def compute(self):
        return NotImplemented


class GenomicIntervalsFull(GenomicIntervals):
    is_stream=False

    def __init__(self, intervals: Interval, chrom_sizes: Dict[str, int], is_stranded=False):
        self._intervals = intervals
        self._geometry = Geometry(chrom_sizes)
        self._chrom_sizes = chrom_sizes
        self._is_stranded = is_stranded

    def __repr__(self):
        return f'Genomic Intervals on {list(self._chrom_sizes)[:5]+["..."]}:\n{self._intervals.astype(Interval)}'

    def get_data(self):
        return self._intervals

    def __len__(self):
        return len(self._intervals)

    def sorted(self):
        return NotImplemented

    def __getitem__(self, idx):
        return self.__class__(self._intervals[idx], self._chrom_sizes, self._is_stranded)

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
            data, self._chrom_sizes, is_stranded=self.is_stranded(),
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
        return self.from_intervals(self._geometry.extend_to_size(self._intervals, size), self._chrom_sizes)

    def merged(self, distance: int = 0) -> GenomicIntervals:
        return self.from_intervals(
            self._geometry.merge_intervals(self._intervals, distance), self._chrom_sizes)

    def get_pileup(self) -> GenomicArray:
        return self._geometry.get_pileup(self._intervals)

    def get_mask(self) -> GenomicArray:
        return self._geometry.get_mask(self._intervals)

    def clip(self) -> 'GenomicIntervalsFull':
        return self.__class__.from_intervals(self._geometry.clip(self._intervals), self._chrom_sizes)

    def __replace__(self, **kwargs):
        return self.__class__(dataclasses.replace(self._intervals, **kwargs), self._chrom_sizes)

    def compute(self):
        return self

    def as_stream(self):
        interval_class = StrandedInterval if self._is_stranded else Interval
        filled = self.genome_context.iter_chromosomes(self._intervals, interval_class)
        # grouped = groupby(self._intervals, 'chromosome')
        # filled = fill_grouped(grouped, self._chrom_sizes.keys(), interval_class)
        return GenomicIntervalsStreamed(
            StreamNode(filled),
            self._chrom_sizes, self._is_stranded)

    def get_sorted_stream(self):
        sorted_intervals = self.sorted()
        return self.from_interval_stream(iter([sorted_intervals]))

    def is_stranded(self):
        return self._is_stranded


class GenomicIntervalsStreamed(GenomicIntervals):
    is_stream = True

    def _get_chrom_size(self, intervals: Interval):
        return self._chrom_sizes[intervals.chromosome]

    def __str__(self):
        return 'GIS:' + str(self._intervals_node)

    def __repr__(self):
        return 'GIS:' + str(self._intervals_node)

    def __init__(self, intervals_node: Node, chrom_sizes: Dict[str, int], is_stranded=False):
        self._chrom_sizes = chrom_sizes
        self._start = ComputationNode(getattr, [intervals_node, 'start'])
        self._stop = ComputationNode(getattr, [intervals_node, 'stop'])
        if is_stranded:
            self._strand = ComputationNode(getattr, [intervals_node, 'strand'])
        self._chromosome = ComputationNode(getattr, [intervals_node, 'chromosome'])
        self._chrom_size_node = StreamNode(iter(self._chrom_sizes.values()))
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
        return self.__class__(ComputationNode(lambda x, i: x[i], [self._intervals_node, item]), self._chrom_sizes)

    def extended_to_size(self, size: int) -> GenomicIntervals:
        return self.__class__(
            ComputationNode(extend_to_size, [self._intervals_node, size, self._chrom_size_node]),
            self._chrom_sizes)
        return self.from_intervals(self._geometry.extend_to_size(self._intervals, size), self._chrom_sizes)

    def merged(self, distance: int = 0) -> GenomicIntervals:
        return self.__class__(ComputationNode(merge_intervals, [self._intervals_node]), self._chrom_sizes)

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
                                self._chrom_sizes)
    
    def get_mask(self) -> GenomicArray:
        return GenomicArrayNode(ComputationNode(get_boolean_mask, [self._intervals_node, self._chrom_size_node]),
                                self._chrom_sizes)

    def clip(self) -> 'GenomicIntervalsFull':
        return self.__class__(ComputationNode(clip, [self._intervals_node, self._chrom_size_node]), self._chrom_sizes)
        return self.__class__.from_intervals(self._geometry.clip(self._intervals), self._chrom_sizes)

    def __replace__(self, **kwargs):
        return self.__class__(ComputationNode(dataclasses.replace, [self._intervals_node], kwargs), self._chrom_sizes)
        return self.__class__(dataclasses.replace(self._intervals, **kwargs), self._chrom_sizes)

    def compute(self):
        chromosome, start, stop = compute(self.chromosome, self.start, self.stop)
        return GenomicIntervalsFull(Interval(chromosome, start, stop), self._chrom_sizes)

    def as_stream(self):
        return self
