from typing import List, Union, Iterable, Tuple, Dict
from abc import abstractclassmethod, abstractmethod, abstractproperty
from npstructures import RunLengthRaggedArray
from .genome_context_base import GenomeContextBase
from ..arithmetics.intervals import GenomicRunLengthArray
from ..datatypes import Interval
import numpy as np
import logging

logger = logging.getLogger(__name__)
GenomeIndex = Union[str, List[str], Interval, Interval.single_entry]


class GenomicData:
    '''
    Base class for genomic data. All genomic data should support indexing on:
    chromosome(s), GenomicIntervals, GenomicLocation and boolean GenomicArrays
    '''

    def __getitem__(self, idx: GenomeIndex):
        if isinstance(idx, str):
            return self.extract_chromsome(idx)
        if (hasattr(idx, 'start') and hasattr(idx, 'stop') and hasattr(idx, 'chromosome')):
            return self.extract_intervals(idx, stranded=hasattr(idx, 'is_stranded') and idx.is_stranded())
        if (hasattr(idx, 'position') and hasattr(idx, 'chromosome')):
            return self.extract_locations(idx, stranded=idx.is_stranded())

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
    def from_stream(cls, stream: Iterable[Tuple[str, GenomicRunLengthArray]], genome_context: GenomeContextBase) -> 'GenomicData':
        return NotImplemented

    @abstractclassmethod
    def from_global_data(cls, global_pileup: GenomicRunLengthArray, genome_context: GenomeContextBase) -> 'GenomicData':
        return NotImplemented

    @abstractmethod
    def to_dict(self) -> Dict[str, np.ndarray]:
        return NotImplemented

    @abstractmethod
    def get_data(self):
        return NotImplemented
