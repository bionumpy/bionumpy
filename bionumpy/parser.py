import gzip
import numpy as np
from npstructures import RaggedArray, RaggedView
from .sequences import Sequences
from .encodings import BaseEncoding, ACTGTwoBitEncoding
from .file_buffers import *

class BufferedNumpyParser:

    def __init__(self, file_obj, buffer_type, chunk_size=1000000):
        self._file_obj = file_obj
        self._chunk_size = chunk_size
        self._is_finished = False
        self._buffer_type = buffer_type

    @classmethod
    def from_filename(cls, filename, *args, **kwargs):
        if any(filename.endswith(suffix) for suffix in ("fasta", "fa", "fasta.gz", "fa.gz")):
            buffer_type = OneLineFastaBuffer
        elif any(filename.lower().endswith(suffix) for suffix in ("fastq", "fq", "fastq.gz", "fq.gz")):
            buffer_type = FastQBuffer
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
