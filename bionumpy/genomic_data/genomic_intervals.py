from abc import ABC, abstractmethod, abstractproperty, abstractclassmethod
from typing import List, Iterable, Tuple, Dict
from .genomic_track import GenomicArray, GenomicArrayNode
from ..datatypes import Interval
from ..arithmetics.intervals import get_pileup, merge_intervals, extend_to_size, clip, get_boolean_mask
from ..streams import groupby
from ..computation_graph import StreamNode, Node, ComputationNode, compute
from .geometry import Geometry
import dataclasses


class GenomicIntervals(ABC):
    @abstractproperty
    def start(self):
        return NotImplemented

    @abstractproperty
    def stop(self):
        return NotImplemented

    @abstractproperty
    def chromosome(self):
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
    def from_intervals(cls, intervals: Interval, chrom_sizes: Dict[str, int]):
        """Create genomic intervals from interval entries and genome info

        Parameters
        ----------
        intervals : Interval
        chrom_sizes : Dict[str, int]
        """
        if isinstance(intervals, Interval): #TODO check is node
            return GenomicIntervalsFull(intervals, chrom_sizes)
        else:
            return cls.from_interval_stream(intervals, chrom_sizes)

    @classmethod
    def from_interval_stream(cls, interval_stream: Iterable[Interval], chrom_sizes: Dict[str, int]):
        """Create streamed genomic intervals from a stream of intervals and genome info

        Parameters
        ----------
        interval_stream : Iterable[Interval]
        chrom_sizes : Dict[str, int]
        """
        
        interval_stream = groupby(interval_stream, 'chromosome')
        interval_stream = StreamNode(pair[1] for pair in interval_stream)
        return GenomicIntervalsStreamed(interval_stream, chrom_sizes)

    @abstractmethod
    def clip(self) -> 'GenomicIntervals':
        return NotImplemented

    def compute(self):
        return NotImplemented


class GenomicIntervalsFull(GenomicIntervals):
    def __init__(self, intervals: Interval, chrom_sizes: Dict[str, int]):
        self._intervals = intervals
        self._geometry = Geometry(chrom_sizes)
        self._chrom_sizes = chrom_sizes

    def __repr__(self):
        return f'Genomic Intervals on {list(self._chrom_sizes)[:5]+["..."]}:\n{self._intervals.astype(Interval)}'

    def get_data(self):
        return self._intervals

    def __len__(self):
        return len(self._intervals)

    def sorted(self):
        return NotImplemented

    def __getitem__(self, idx):
        return self.__class__(self._intervals[idx], self._chrom_sizes)

    @property
    def start(self):
        return self._intervals.start

    @property
    def stop(self):
        return self._intervals.stop

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


class GenomicIntervalsStreamed:
    def _get_chrom_size(self, intervals: Interval):
        return self._chrom_sizes[intervals.chromosome]

    def __str__(self):
        return 'GIS:' + str(self._intervals_node)

    def __init__(self, intervals_node: Node, chrom_sizes: Dict[str, int]):
        self._chrom_sizes = chrom_sizes
        self._start = ComputationNode(getattr, [intervals_node, 'start'])
        self._stop = ComputationNode(getattr, [intervals_node, 'stop'])
        # self._strand = ComputationNode(getattr, [intervals_node, 'strand'])
        self._chromosome = ComputationNode(getattr, [intervals_node, 'chromosome'])
        self._chrom_size_node = StreamNode(iter(self._chrom_sizes.values()))
        self._intervals_node = intervals_node

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
        return GenomicIntervalsFull(Interval(chromosome, start, stop))
