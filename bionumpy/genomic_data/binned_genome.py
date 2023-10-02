import numpy as np
from bionumpy import LocationEntry, bnp_open
from bionumpy.genomic_data.genome import Genome
from bionumpy.genomic_data.genome_context_base import GenomeContextBase


class BinnedGenome:
    def __init__(self, genome_context: GenomeContextBase, bin_size: int = 1000):
        self._genome_context = genome_context
        self._bin_size = bin_size
        chrom_sizes = np.array(list(genome_context.chrom_sizes.values()))
        self._n_bins = (chrom_sizes+bin_size-1)//bin_size
        self._bin_offsets = np.insert(np.cumsum(self._n_bins), 0, 0)
        self._counts = np.zeros(self._bin_offsets[-1], dtype=int)

    @classmethod
    def from_file(cls, filename: str, bin_size: int = 1000):
        genome = Genome.from_file(filename)
        return cls(genome.get_genome_context(), bin_size)

    @property
    def genome_context(self):
        return self._genome_context

    def count(self, entries: LocationEntry, position_field='position'):
        chrom_nrs = self._genome_context.encoding.encode(entries.chromosome).raw()
        offsets = getattr(entries, position_field)//self._bin_size
        bin_nr = self._bin_offsets[chrom_nrs] + offsets
        self._counts += np.bincount(bin_nr, minlength=self._bin_offsets[-1])

    def count_file(self, filename, position_field='position'):
        chunks = bnp_open(filename, 'r').read_chunks()
        for chunk in chunks:
            self.count(chunk, position_field=position_field)

    @property
    def count_dict(self):
        return {chrom: self._counts[self._bin_offsets[i]:self._bin_offsets[i+1]]
                for i, chrom in enumerate(self._genome_context.chrom_sizes)}

    def __getitem__(self, chromosome):
        i = self._genome_context.encoding.encode(chromosome).raw()
        return self._counts[self._bin_offsets[i]:self._bin_offsets[i+1]]