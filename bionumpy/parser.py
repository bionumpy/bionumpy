import time
import logging
import sys
from shared_memory_wrapper import SingleSharedArray, to_shared_memory, from_shared_memory
import gzip
import numpy as np

from .sequences import Sequences
from .encodings import BaseEncoding, ACTGTwoBitEncoding
NEWLINE = 10

def get_mask_from_intervals(intervals, size):
    """ intervals = (starts, ends) """
    mask_changes = np.zeros(size+1, dtype="bool")
    mask_changes[intervals[0]]=True
    mask_changes[intervals[1]]^=True
    mask = np.logical_xor.accumulate(mask_changes)
    return mask[:-1]

class FileBuffer:
    _buffer_divisor = 1
    COMMENT = 0
    def __init__(self, data, new_lines):
        self._data = data
        self._new_lines = new_lines
        self._is_validated = False
        self.size = self._data.size

    @classmethod
    def from_raw_buffer(cls, chunk):
        raise NotImplemented


    def validate_if_not(self):
        if not self._is_validated:
            self._validate()
    
    def _move_intervals_to_contigous_array(self, starts, ends):
        """
        Create a mask for where the sequences are and move sequences to continous array
        """
        mask = get_mask_from_intervals((starts, ends), self._data.size)
        removed_areas = np.cumsum(starts-np.insert(ends[:-1], 0, 0))
        new_intervals = (starts-removed_areas, ends-removed_areas)
        m = np.sum(mask)
        d = m % self._buffer_divisor
        seq = np.empty(m-d+self._buffer_divisor, dtype=self._data.dtype)
        seq[:m] = self._data[mask]
        return seq, new_intervals

    def _move_intervals_to_2d_array(self, starts, ends, fill_value=0):
        n_intervals = starts.size
        n_chars = ends-starts
        max_chars = np.max(n_chars)
        array = np.full((n_intervals, max_chars), fill_value, dtype=np.uint8)
        from_mask = get_mask_from_intervals((starts, ends), self._data.size)
        to_ends = max_chars*np.arange(1, n_intervals+1)
        to_starts = to_ends-n_chars
        to_mask = get_mask_from_intervals((to_starts, to_ends), array.size).reshape(array.shape)
        array[to_mask] = self._data[from_mask]
        return array


class OneLineBuffer(FileBuffer):
    n_lines_per_entry = 2
    _encoding = BaseEncoding

    @classmethod
    def from_raw_buffer(cls, chunk):
        new_lines = np.flatnonzero(chunk==NEWLINE)
        n_lines = new_lines.size
        assert n_lines >= cls.n_lines_per_entry, "No complete entry in buffer"
        new_lines = new_lines[:n_lines-(n_lines%cls.n_lines_per_entry)]
        return cls(chunk[:new_lines[-1]+1], new_lines)

    def get_sequences(self):
        self.validate_if_not()
        sequence_starts = self._new_lines[::self.n_lines_per_entry]+1
        sequence_ends = self._new_lines[1::self.n_lines_per_entry]
        
        seq, new_intervals = self._move_intervals_to_contigous_array(sequence_starts, sequence_ends)
        seq = self._encoding.from_bytes(seq)
        print(seq, new_intervals)
        lens = new_intervals[1]-new_intervals[0]
        return Sequences(seq, lens, encoding=self._encoding)

    def _validate(self):
        n_lines = self._new_lines.size
        assert n_lines % self.n_lines_per_entry == 0, "Wrong number of lines in buffer"
        header_idxs = self._new_lines[self.n_lines_per_entry-1:-1:self.n_lines_per_entry]+1
        assert np.all(self._data[header_idxs]==self.HEADER)
        self._is_validated = True

class OneLineFastaBuffer(OneLineBuffer):
    HEADER= 62
    n_lines_per_entry = 2

class FastQBuffer(OneLineBuffer):
    HEADER= 64
    n_lines_per_entry = 4
    _encoding = BaseEncoding

class OneLineFastaBuffer2Bit(OneLineFastaBuffer):
    _encoding = ACTGTwoBitEncoding
    _buffer_divisor = 32

class FastQBuffer2Bit(FastQBuffer):
    _encoding = ACTGTwoBitEncoding
    _buffer_divisor = 32

class BufferedNumpyParser:

    def __init__(self, file_obj, buffer_type, chunk_size=1000000):
        self._file_obj = file_obj
        self._chunk_size = chunk_size
        self._is_finished = False
        self._buffer_type = buffer_type

    @classmethod
    def from_filename(cls, filename, *args, **kwargs):
        if any(filename.endswith(suffix) for suffix in ("fasta", "fa", "fasta.gz", "fa.gz")):
            buffer_type = OneLineFastaBuffer2Bit
        elif any(filename.lower().endswith(suffix) for suffix in ("fastq", "fq", "fastq.gz", "fq.gz")):
            buffer_type = FastQBuffer2Bit
        else:
            raise NotImplemented
        if filename.endswith(".gz"):
            file_obj = gzip.open(filename, "rb")
        else:
            file_obj = open(filename, "rb")
        return cls(file_obj, buffer_type, *args, **kwargs)

    def get_chunk(self):
        a, bytes_read = self.read_raw_chunk()
        self._is_finished = bytes_read < self._chunk_size
        if bytes_read == 0:
            return None
        
        # Ensure that the last entry ends with newline. Makes logic easier later
        if self._is_finished  and a[bytes_read-1] != NEWLINE:
            a[bytes_read] = NEWLINE
            bytes_read += 1
        return a[:bytes_read]

    def read_raw_chunk(self):
        array = np.empty(self._chunk_size, dtype="uint8")
        bytes_read = self._file_obj.readinto(array)
        return array, bytes_read

    def remove_initial_comments(self):
        if self._buffer_type.COMMENT == 0:
            return 
        for line in self._file_obj:
            if line[0] != self._buffer_type.COMMENT:
                self._file_obj.seek(-len(line), 1)
                break

    def get_chunks(self):
        self.remove_initial_comments()
        chunk = self.get_chunk()
        buff = self._buffer_type.from_raw_buffer(chunk)
        while not self._is_finished:
            buff = self._buffer_type.from_raw_buffer(chunk)
            self._file_obj.seek(buff.size-self._chunk_size, 1)
            yield buff
            chunk = self.get_chunk()
        if chunk is not None and chunk.size:
            yield self._buffer_type.from_raw_buffer(chunk)

class TextParser:
    _encoding = BaseEncoding

    def __init__(self, filename, chunk_size=1000000):
        self._file_obj = self._get_file_obj(filename)
        self._chunk_size = chunk_size
        self._is_finished = False
        self._last_header_idx = -1

    def _get_file_obj(self, filename):
        if filename.endswith(".gz"):
            return gzip.open(filename, "rb")
        return open(filename, "rb")

    def _cut_array(self, array, last_newline):
        self._last_header_idx = last_newline+1
        return array[:last_newline]

    def read_raw_chunk(self):
        array = np.empty(self._chunk_size, dtype="uint8")
        bytes_read = self._file_obj.readinto(array)
        return array, bytes_read

    def parse(self, as_shared_memory_object=False):
        while not self._is_finished:
            t = time.perf_counter()
            chunk = self.parse_chunk()
            logging.info("Time to parse chunk: %.4f" % (time.perf_counter()-t))
            if chunk is not None:
                if as_shared_memory_object:
                    shared_memory_name = str(np.random.randint(0, 10e15))
                    to_shared_memory(chunk, shared_memory_name)
                    yield shared_memory_name
                else:
                    yield chunk

    def _handle_chunk(self, chunk):
        return self.get_sequences(chunk)

    def get_raw_chunks(self, as_shared_memory=False):
        raw_chunk = self.read_raw_chunk()
        while not self._is_finished and raw_chunk:
            if as_shared_memory_object:
                shared_memory_name = str(np.random.randint(0, 10e15))
                to_shared_memory(chunk, shared_memory_name)
                yield shared_memory_name
            else:
                yield chunk

    def cut_chunk(self, chunk):
        new_lines = np.flatnonzero(chunk == NEW_LINE)
        return self._cut_array(chunk, new_lines[-1])
        

    def parse_chunk(self):
        a, bytes_read = self.read_raw_chunk()
        self._is_finished = bytes_read < self._chunk_size
        if bytes_read == 0:
            return None
        if self._is_finished  and a[bytes_read-1] != NEWLINE:
            a[bytes_read]  = NEWLINE
        sequences = self._handle_chunk(a[:bytes_read])
        if not self._is_finished:
            self._file_obj.seek(self._last_header_idx-self._chunk_size, 1)
        return sequences

class OneLineParser(TextParser):
    """
    Base class for formats where the sequnences
    are constrained to sinlge line
    """
    _buffer_divisor=1

    def cut_chunk(self, chunk):
        new_lines = np.flatnonzero(chunk == NEW_LINE)
        return self._cut_array(chunk, new_lines[-1])

    def get_sequences(self, raw_chunk):
        new_lines = np.flatnonzero(raw_chunk==NEWLINE)
        chunk, new_lines = self._validate_chunk(raw_chunk, new_lines)
        sequence_starts = new_lines[::self.n_lines_per_entry]+1
        sequence_ends = new_lines[1::self.n_lines_per_entry]
        return self._mask_and_move_sequences(chunk, sequence_starts, sequence_ends)

    def _mask_and_move_sequences(self, array, sequence_starts, sequence_ends):
        """
        Create a mask for where the sequences are and move sequences to continous array
        """
        mask = get_mask_from_intervals((sequence_starts, sequence_ends), array.size)
        removed_areas = np.cumsum(sequence_starts-np.insert(sequence_ends[:-1], 0, 0))
        new_intervals = (sequence_starts-removed_areas, sequence_ends-removed_areas)
        m = np.sum(mask)
        d = m%self._buffer_divisor
        seq = np.empty(m-d+self._buffer_divisor, dtype=array.dtype)
        seq[:m] = array[mask]
        seq = self._encoding.from_bytes(seq)
        return Sequences(seq, new_intervals, encoding=self._encoding)

    def _validate_chunk(self, chunk, new_lines):
        assert chunk[0] == self.HEADER, "Chunk does not start with header"
        n_lines = new_lines.size
        assert n_lines >=self.n_lines_per_entry, "No complete entry in buffer"
        new_lines = new_lines[:n_lines-(n_lines%self.n_lines_per_entry)]
        header_idxs = new_lines[self.n_lines_per_entry-1:-1:self.n_lines_per_entry]+1
        assert np.all(chunk[header_idxs]==self.HEADER)
        return self._cut_array(chunk, new_lines[-1]), new_lines

class FastqParser(OneLineParser):
    HEADER = 64
    n_lines_per_entry = 4


class FastqParser2Bit(FastqParser):
    _buffer_divisor = 32
    _encoding = ACTGTwoBitEncoding


class OneLineFastaParser(OneLineParser):
    HEADER= 62
    n_lines_per_entry = 2


class LooseOneLineFastaParser(OneLineParser):
    HEADER = 62
    n_lines_per_entry = 2

    def _validate_chunk(self, chunk, new_lines):
        assert chunk[0]==self.HEADER

    def get_sequences(self, array):
        new_lines = np.flatnonzero(array==NEWLINE)
        # Find which lines are sequences, and cut the array so that it ends with a complete sequence
        is_sequence_line = array[new_lines[:-1]+1] != self.HEADER
        # assert np.any(is_sequence_line), (array, self._is_finished)
        idx_last_sequence_line = np.flatnonzero(is_sequence_line)[-1]
        new_lines = new_lines[:idx_last_sequence_line+2]
        array = self._cut_array(array, new_lines[-1])
        is_sequence_line = is_sequence_line[:idx_last_sequence_line+1]
        sequence_starts = new_lines[:-1][is_sequence_line]+1
        sequence_ends = new_lines[1:][is_sequence_line]
        return self._mask_and_move_sequences(array, sequence_starts, sequence_ends)

    def _mask_and_move_sequences__(self, array, sequence_starts, sequence_ends):
        """
        Create a mask for where the sequences are and move sequences to continous array
        """
        mask = get_mask_from_intervals((sequence_starts, sequence_ends), array.size)
        removed_areas = np.cumsum(sequence_starts-np.insert(sequence_ends[:-1], 0, 0))
        new_intervals = (sequence_starts-removed_areas, sequence_ends-removed_areas)
        return Sequences(array[mask], new_intervals[0], new_intervals[1])


class OneLineFastaParser2bit(OneLineFastaParser):
    _buffer_divisor = 32
    _encoding = ACTGTwoBitEncoding
    def _mask_and_move_sequences__(self, array, sequence_starts, sequence_ends):
        """
        Create a mask for where the sequences are and move sequences to continous array
        """
        mask = get_mask_from_intervals((sequence_starts, sequence_ends), array.size)
        removed_areas = np.cumsum(sequence_starts-np.insert(sequence_ends[:-1], 0, 0))
        new_intervals = (sequence_starts-removed_areas, sequence_ends-removed_areas)
        m = np.sum(mask)
        d = m%32
        seq = np.empty(m-d+32, dtype=array.dtype)
        seq[:m] = array[mask]
        return Sequences(seq, new_intervals[0], new_intervals[1])


class FastaParser(TextParser):
    def _get_sequence_offsets(self, new_lines, is_header):
        header_starts = new_lines[:-1][is_header[:-1]]
        header_ends = new_lines[1:][is_header[:-1]]
        header_lens = header_ends-header_starts
        masked_count = np.ones(new_lines.size, dtype="int")
        masked_count[1:][is_header[:-1]] = header_lens
        masked_count[0] = new_lines[0]+1
        cum_masked_count = np.cumsum(masked_count)
        indexes_after_masking = new_lines+1-cum_masked_count
        starts = indexes_after_masking[1:][is_header[:-1]]
        if not self._is_finished:
            starts = starts[:-1]
        return starts

    def get_sequences(self, array):
        assert array[0] == HEADER
        all_new_lines = np.flatnonzero(array==NEWLINE)
        new_lines = all_new_lines if all_new_lines[-1] != array.size-1 else all_new_lines[:-1]
        is_header = array[new_lines+1] == HEADER
        mask_changes = np.zeros(array.size, dtype="bool")
        mask_changes[new_lines+1] = np.where(is_header, 1, 0)
        mask_changes[new_lines[1:]] = np.where(is_header[:-1], 1, 0)
        mask_changes[0] = 1
        mask_changes[new_lines[0]] = 1
        mask = np.logical_xor.accumulate(mask_changes)
        mask[all_new_lines] = 1
        beginnings = self._get_sequence_offsets(new_lines, is_header)
        intervals = (np.insert(beginnings, 0, 0), np.append(beginnings, np.count_nonzero(~mask)))
        if self._is_finished:
            return Sequences(array[~mask], intervals)
        assert np.any(is_header), ("Fasta entry to big for chunk size", to_text(array), chunk_size)
        self._last_header_idx = new_lines[is_header][-1]+1
        intervals = (np.insert(beginnings, 0, 0), np.append(beginnings, np.count_nonzero(~mask[:self._last_header_idx])))
        return Sequences(array[:self._last_header_idx][~mask[:self._last_header_idx]], intervals)
