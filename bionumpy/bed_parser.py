from .parser import FileBuffer, NEWLINE, get_mask_from_intervals
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

    def _validate(self):
        chunk = self._data
        delimiters = self._delimiters[1:]
        n_delimiters_per_line = next(i for i, d in enumerate(delimiters) if chunk[d] == NEWLINE) + 1
        self._n_cols = n_delimiters_per_line
        last_new_line = next(i for i, d in enumerate(delimiters[::-1]) if chunk[d] == NEWLINE)
        delimiters = delimiters[:delimiters.size-last_new_line]
        assert delimiters.size % n_delimiters_per_line == 0, "irregular number of delimiters per line"
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
        assert data.shape[-1] == 2
        
        self.start = self.data[..., 0]
        self.end = self.data[..., 1]

    def __repr__(self):
        return repr(self.data)

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
        idx = np.searchsorted(self.starts, position)
        return position < self.ends[idx]

class BedBuffer(DelimitedBuffer):
    def get_intervals(self):
        self.validate_if_not()
        data =  self.get_integers(cols=[1, 2])
        return Interval(data)

class ChromosomeProvider:
    @staticmethod
    def get_chrom_name(char_array):
        return "".join(chr(c) for c in char_array).replace("\x00", "")
        
    @staticmethod
    def _is_same_chromosome(chrom_1, chrom_2):
        return FullBedFile.get_chrom_name(chrom_1) == FullBedFile.get_chrom_name(chrom_2)

    @staticmethod
    def _get_chromosome_changes(chromosomes):
        return np.flatnonzero(
            np.all(chromosomes[1:] != chromosomes[:-1], axis=-1))+1

class ChromosomeDictProvider:
    pass

class ChromosomeStreamProvider:
    def __init__(self, file_buffers):
        self._buffers = file_buffers

    def __iter__(self):
        cur_data = []
        last_chromosome = np.zeros(3, dtype=np.uint8)
        for file_buffer in self._buffers:
            chromosomes = file_buffer.get_chromosomes()
            if not len(chromosomes):
                break
            if not cls._is_same_chromosome(last_chromosome, chromosomes[0]):
                yield (cls.get_chrom_name(last_chromosome), np.concatenate(cur_data))
                last_chromosome = chromosomes[0]
                cur_data = []
            data = file_buffer.get_data()
            chromosome_changes = cls._get_chromosome_changes(chromosomes)
            if len(chromosome_changes)==0:
                cur_data.append(data)
                
            cur_intervals.append(data[:chromosome_changes[0]])
            yield np.get_chrom_name(last_chromosome), np.concatenate(cur_intervals)
            for start, end in zip(chromosome_changes[:-1], chromosome_changes[1:]):
                yield np.get_chrom_name(chromosomes[start]), data[start:end]
            last_chromosome = chromosomes[-1]
            cur_data.append(data[chromosome_changes[-1]:])
        yield cls.get_chrom_name(last_chromosome). np.concatenate(cur_data)

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
#     def __init__(self, chromosome, position, ref_seq, alt_seq):
#         self._chromosome = chromosome
#         self._position = position
#         self._ref_seq = ref_seq
#         self._alt_seq = alt_seq

class VCFBuffer(DelimitedBuffer):
    def get_snps(self):
        self.validate_if_not()
        chromosomes = self.get_text(0)
        position = self.get_integers(1)
        from_seq = self.get_text(3)
        to_seq = self.get_text(4)
        return SNP(chromosomes, position, from_seq, to_seq)

def apply_to_chromosomes(func):
    def new_func(*args, **kwargs):
        assert all(isinstance(arg, ChromosomeProvider) for arg in chain(args, kwargs.values()))
        args_streams = [isinstance(a, ChromosomeStreamProvider) for a in args]
        kwargs_streams = {key: isinstance(value, ChromosomeStreamProvider) for key, value in kwargs.items()}
        n_streams = sum(chain(args_streams, kwarg_streams))
        assert sum(chain(args_streams, kwarg_streams)) <= 1
        if n_streams == 0:
            all_chroms = {key for provider in chain(args, kwargs.values()) for   key in proivder.keys()}
            sorted_chroms = sorted(all_chroms)
            for chromsome in sorted_chroms:
                new_args = [arg[chrom] for arg in args]

                new_kwargs = {key: val[chrom] for kwy, val in kwargs.items()}
                yield (chrom, func(*new_args, **new_kwargs)
        else:
            stream = next(chain(
                (arg for is_stream, arg in zip(args_streams, args) if is_stream),
                (val for key, val in kwargs.items() if kwarg_streams[key])))

            for chromosome, data in stream:
                new_args = [data if is_stream else arg[chromosome] for is_stream, arg in zip(arg_streams, args)]
                nww_kwargs = {key: data if kwargs_streams[key] else value[chromosome]
                              for key, value in kwargs.items()}
                yield chromosome, func(*new_args, **new_kwargs)
    return new_func
