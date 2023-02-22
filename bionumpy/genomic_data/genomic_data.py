from typing import List, Union, Iterable, Tuple, Dict, Any
from abc import abstractclassmethod, abstractmethod, abstractproperty, ABC
from npstructures import RunLengthRaggedArray
from ..arithmetics.intervals import GenomicRunLengthArray
from .global_offset import GlobalOffset
from ..datatypes import  Interval
from ..streams import groupby
import numpy as np
import logging

logger = logging.getLogger(__name__)
GenomeIndex = Union[str, List[str], Interval, Interval.single_entry]


def fill_grouped(grouped: Iterable[Tuple[str, Any]], real_order: Iterable[str], dataclass: type, allowed: set = None):
    real_order = iter(real_order)
    next_real = next(real_order, None)
    for name, group in grouped:
        assert next_real is not None
        while name != next_real:
            print('yielding empty', name, dataclass)
            yield dataclass.empty()
            next_real = next(real_order, None)
        print('yielding full', name, dataclass)
        yield group
        next_real = next(real_order, None)
    # next_real = next(real_order, None)
    while next_real is not None:
        print('yielding last', next_real)
        yield next_real
        next_real = next(real_order, None)


class GenomeError(Exception):
    pass


class GenomeContext:
    def __init__(self, chrom_size_dict: Dict[str, int], ignored=None):
        self._chrom_size_dict = chrom_size_dict
        self._ignored = ignored
        if ignored is None:
            self._ignored = set([])
        self._included = set(chrom_size_dict.keys())

    @property
    def chrom_sizes(self):
        return self._chrom_size_dict

    @classmethod
    def from_dict(cls, chrom_size_dict):
        f = lambda key: '_' not in key
        return cls({key: value for key, value in chrom_size_dict.items() if f(key)},
                   {key for key in chrom_size_dict if not f(key)})

    def chromosome_order(self):
        return (key for key in self._chrom_size_dict if '_' not in key)

    def is_compatible(self, other):
        return self._chrom_size_dict == other.chrom_size_dict

    def _included_groups(self, grouped):
        for name, group in grouped:
            if name in self._ignored:
                continue
            if name not in self._included:
                raise GenomeError(f'{name} not included in genome: {set(self._chrom_size_dict.keys())}')
            yield name, group

    def iter_chromosomes(self, data, dataclass, group_field='chromosome'):
        real_order = self.chromosome_order()# self._chrom_size_dict.keys()
        grouped = groupby(data, group_field)
        grouped = self._included_groups(grouped)
        next_name, next_group = next(grouped, (None, None))
        seen = []
        seen_group = []
        for name in real_order:
            if name == next_name:
                seen_group.append(next_name)
                logger.debug(f'Yielding data for {name}')
                yield next_group
                next_name, next_group = next(grouped, (None, None))
                if next_name in seen:
                    raise GenomeError(f'Sort order discrepancy ({next_name}): Genome so far {seen}, data: {seen_group}')
            else:
                logger.debug(f'Yielding empty data for {name}')
                yield dataclass.empty()
            seen.append(name)
        if next(grouped, None) is not None:
            raise GenomeError()
        # 
        # real_order = iter(real_order)
        # 
        # next_real = next(real_order, None)
        # for name, group in self._included_groups(grouped):
        #     assert next_real is not None
        #     while name != next_real:
        #         print('yielding empty', name, dataclass)
        #         yield dataclass.empty()
        #         next_real = next(real_order, None)
        #     print('yielding full', name, dataclass)
        #     yield group
        #     next_real = next(real_order, None)
        # # next_real = next(real_order, None)
        # while next_real is not None:
        #     print('yielding last', next_real)
        #     yield next_real
        #     next_real = next(real_order, None)


class GenomicData:
    def __getitem__(self, idx: GenomeIndex):
        if isinstance(idx, str):
            return self.extract_chromsome(idx)
        if (hasattr(idx, 'start') and hasattr(idx, 'stop') and hasattr(idx, 'chromosome') and hasattr(idx, 'is_stranded')):
            return self.extract_intervals(idx, stranded=idx.is_stranded())
        if isinstance(idx, list):
            if len(idx) == 0:
                return self.empty()
            if isinstance(idx[0], str):
                return self.extract_chromosome(idx)
        if isinstance(idx, GenomicData) and idx.dtype == bool:
            return self._index_boolean(idx)
        raise ValueError(f'{type(idx)} object not valid as index for GenomicData: {idx}')

    @abstractproperty
    def genome_context(self):
        return NotImplemented

    @abstractmethod
    def _index_boolean(self, chromosome: Union[str, List[str]]) -> 'GenomicData':
        return NotImplemented

    def dtype(self):
        return None

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
