from .parser import FileBuffer, NEWLINE
from .chromosome_provider import *
from dataclasses import dataclass
import numpy as np

class DelimitedBuffer(FileBuffer):
    DELIMITER = ord("\t")
    COMMENT = ord("#")

    def __init__(self, data, new_lines):
        super().__init__(data, new_lines)
        self._delimiters = np.concatenate(([0],
            np.flatnonzero(self._data == self.DELIMITER),
            self._new_lines))
        self._delimiters.sort(kind="mergesort")

    @classmethod
    def from_raw_buffer(cls, chunk):
        new_lines = np.flatnonzero(chunk==NEWLINE)
        return cls(chunk[:new_lines[-1]+1], new_lines)

    def get_integers(self, cols):
        cols = np.asanyarray(cols)
        integer_starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, cols]+1
        integer_ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, cols]
        integers = self._extract_integers(integer_starts.ravel(), integer_ends.ravel())
        return integers.reshape(-1, cols.size)

    def _extract_integers(self, integer_starts, integer_ends):
        digit_chars = self._move_intervals_to_2d_array(integer_starts, integer_ends, DigitEncoding.MIN_CODE)
        n_digits = digit_chars.shape[-1]
        powers = np.uint32(10)**np.arange(n_digits)[::-1]
        return DigitEncoding.from_bytes(digit_chars) @ powers

    def get_text(self, col):
        self.validate_if_not()
        # delimiters = self._delimiters.reshape(-1, self._n_cols)
        starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, col]+1
        ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, col]
        return self._move_intervals_to_2d_array(starts, ends)

    def get_text_range(self, col, start=0, end=None):
        self.validate_if_not()
        # delimiters = self._delimiters.reshape(-1, self._n_cols)
        starts = self._delimiters[:-1].reshape(-1, self._n_cols)[:, col]+1+start
        if end is not None:
            ends = starts+end
        else:
            ends = self._delimiters[1:].reshape(-1, self._n_cols)[:, col]
        return self._move_intervals_to_2d_array(starts.ravel(), ends.ravel())

    def _validate(self):
        chunk = self._data
        delimiters = self._delimiters[1:]
        n_delimiters_per_line = next(i for i, d in enumerate(delimiters) if chunk[d] == NEWLINE) + 1
        self._n_cols = n_delimiters_per_line
        last_new_line = next(i for i, d in enumerate(delimiters[::-1]) if chunk[d] == NEWLINE)
        delimiters = delimiters[:delimiters.size-last_new_line]
        assert delimiters.size % n_delimiters_per_line == 0, f"irregular number of delimiters per line ({delimiters.size}, {n_delimiters_per_line})"
        delimiters = delimiters.reshape(-1, n_delimiters_per_line)
        assert np.all(chunk[delimiters[:, -1]] == NEWLINE)
        self._validated = True
        
class StrandEncoding:
    MIN_CODE = ord("+")
    @classmethod
    def from_bytes(cls, bytes_array):
        return (bytes_array & np.uint8(2)) >> np.uint8(1)

    @classmethod
    def to_bytes(cls, strands):
        return 2*strands + cls.MIN_CODE

class DigitEncoding:
    MIN_CODE = ord("0")
    @classmethod
    def from_bytes(cls, bytes_array):
        return bytes_array-cls.MIN_CODE

    @classmethod
    def to_bytes(cls, digits):
        return digits+cls.MIN_CODE

class Interval:
    def __init__(self, data):
        self.data = np.asanyarray(data)
        assert self.data.shape[-1] == 2
        
        self.start = self.data[..., 0]
        self.end = self.data[..., 1]

    def __repr__(self):
        return repr(self.data)

    def __eq__(self, other):
        return np.all(self.data==other.data)

    def in_interval(self, position):
        return (self.start <= position) & (position < self.end)

class SortedIntervals:
    def __init__(self, data):
        self.data = np.asanyarray(data)
        assert data.shape[-1] == 2
        assert len(data.shape) == 2
        
        self.starts = self.data[..., 0]
        self.ends = self.data[..., 1]

    def in_intervals(self, position):
        idx = np.minimum(np.searchsorted(self.starts, position, side="left"), self.starts.size-1)
        return (position >=self.starts[idx]) & (position < self.ends[idx])

    @classmethod
    def concatenate(cls, elements):
        return cls(np.vstack([element.data for element in elements]))

class BedBuffer(DelimitedBuffer):
    data_class=SortedIntervals
    def get_intervals(self):
        self.validate_if_not()
        data = self.get_data()
        return Interval(data)

    def get_data(self):
        self.validate_if_not()
        return self.get_integers(cols=[1, 2])

    def get_chromosomes(self):
        self.validate_if_not()
        return self.get_text(col=0)

class FullBedFile(ChromosomeDictProvider):
    def __init__(self, chrom_dict, all_intervals):
        self._chrom_dict = chrom_dict
        self._all_intervals = all_intervals
        
    def __getitem__(self, chrom):
        start, end = chrom_dict[chrom]
        return all_intervals[start:end]

    @classmethod
    def from_bed_buffer_stream(cls, bed_buffer_stream):
        cur_offset = 0
        all_intervals = []
        chrom_starts = {}
        chrom_ends = {}
        for bed_buffer in bed_buffer_stream:
            chromosomes = cur_buffer.get_chromosomes()
            if not len(chromosomes):
                break
            if not cls._is_same_chromosome(last_chromosome, chromosomes[0]):
                chrom_starts[cls.get_chrom_name(chromosomes[0])] == cur_offset
            chromosome_changes = cls._get_chromosome_changes(chromosomes)
            chrom_starts.update((cls.get_chrom_name(chromosomes[i]), i+cur_offset)
                                for i in chromosome_changes)
            all_intervals.append(cur_buffer.get_intervals())
            last_chromosome = chromosomes[-1]
            cur_offset += len(chromosomes)
        chrom_ends[cls.get_chrom_name(last_chromosome)] == cur_offset
        interval_dict = {chrom: (chrom_starts[chrom], chrom_ends[chrom]) for chrom in chrom_starts.keys()}
        return cls(interval_dict, np.concatenate(intervals))

@dataclass
class SNP:
    chromosome: np.array
    position: np.array
    ref_seq: np.array
    alt_seq: np.array

    def __iter__(self):
        return (SNP(*comb) for comb in zip(self.chromosome, self.position, self.ref_seq, self.alt_seq))

    def __getitem__(self, idx):
        return SNP(self.chromosome[idx], self.position[idx], self.ref_seq[idx], self.alt_seq[idx])

#     def __init__(self, chromosome, position, ref_seq, alt_seq):
#         self._chromosome = chromosome
#         self._position = position
#         self._ref_seq = ref_seq
#         self._alt_seq = alt_seq

class VCFBuffer(DelimitedBuffer):
    def get_snps(self):
        self.validate_if_not()
        chromosomes = self.get_text(0)
        position = self.get_integers(1).ravel()-1
        from_seq = self.get_text(3).ravel()
        to_seq = self.get_text(4).ravel()
        return SNP(chromosomes, position, from_seq, to_seq)

    def get_data(self):
        self.validate_if_not()
        chromosomes = self.get_text(0)
        position = self.get_integers(1)
        from_seq = self.get_text(3)
        to_seq = self.get_text(4)
        return SNP(chromosomes, position, from_seq.ravel(), to_seq.ravel())

class GenotypeEncoding:
    @classmethod
    def from_bytes(cls, bytes_array):
        assert bytes_array.shape[-1]==3
        return (bytes_array[..., 0]==ord("1"))+(bytes_array[..., 2]==ord("1")).astype("int")

class VCFMatrixBuffer(VCFBuffer):
    def get_entries(self):
        self.validate_if_not()
        chromosomes = self.get_text(0)
        position = self.get_integers(1).ravel()-1
        from_seq = self.get_text(3).ravel()
        to_seq = self.get_text(4).ravel()
        n_samples = self._n_cols-9
        genotypes = self.get_text_range(np.arange(9, self._n_cols), end=3)
        return SNP(chromosomes, position, from_seq, to_seq), GenotypeEncoding.from_bytes(genotypes.reshape(-1, n_samples, 3))

    def get_data(self):
        self.validate_if_not()
        chromosomes = self.get_text(0)
        position = self.get_integers(1)
        from_seq = self.get_text(3)
        to_seq = self.get_text(4)
        return SNP(chromosomes, position, from_seq.ravel(), to_seq.ravel())
