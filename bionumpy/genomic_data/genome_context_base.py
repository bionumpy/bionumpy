from abc import ABC, abstractproperty, abstractclassmethod, abstractmethod


class GenomeContextBase(ABC):
    @abstractproperty
    def chrom_sizes(self):
        return NotImplemented

    # @abstractproperty
    # def geometry(self):
    #     return NotImplemented

    @abstractproperty
    def global_offset(self):
        return NotImplemented

    @abstractmethod
    def is_included(self, chromosomes):
        return NotImplemented

    @abstractclassmethod
    def from_dict(cls, chrom_size_dict):
        return NotImplemented

    @abstractmethod
    def chromosome_order(self):
        return NotImplemented

    @abstractmethod
    def is_compatible(self, other):
        return NotImplemented

    @abstractmethod
    def iter_chromosomes(self, data, dataclass, group_field='chromosome'):
        return NotImplemented
