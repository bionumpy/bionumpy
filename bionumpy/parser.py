import gzip
import numpy as np
from pathlib import PurePath
from npstructures import RaggedArray, RaggedView
from .sequences import Sequences
from .encodings import BaseEncoding, ACTGTwoBitEncoding
from .file_buffers import *
from .delimited_buffers import *


class BufferedNumpyParser:
    buffer_types = {".vcf": VCFBuffer,
                    ".bed": BedBuffer,
                    ".fasta": TwoLineFastaBuffer,
                    ".fa": TwoLineFastaBuffer,
                    ".fastq": FastQBuffer,
                    ".fq": FastQBuffer}


    def __init__(self, file_obj, buffer_type, chunk_size=1000000):
        self._file_obj = file_obj
        self._chunk_size = chunk_size
        self._is_finished = False
        self._buffer_type = buffer_type

    @classmethod
    def from_filename(cls, filename, *args, **kwargs):
        path = PurePath(filename)
        suffixes = path.suffixes

        open_func = open
        if suffixes[-1] == ".gz":
            open_func = gzip.open
            suffixes = suffixes[:-1]
        buffer_type = cls.buffer_types[suffixes[-1]]
        return cls(open_func(filename, "rb"),
                   buffer_type,
                   *args, **kwargs)

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

    def __iter__(self):
        return self.get_chunks()

    def get_chunks(self):
        self.remove_initial_comments()
        chunk = self.get_chunk()
        # buff = self._buffer_type.from_raw_buffer(chunk)
        while not self._is_finished:
            buff = self._buffer_type.from_raw_buffer(chunk)
            self._file_obj.seek(buff.size-self._chunk_size, 1)
            yield buff
            chunk = self.get_chunk()
        if chunk is not None and chunk.size:
            yield self._buffer_type.from_raw_buffer(chunk)
