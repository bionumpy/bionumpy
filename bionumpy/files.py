from pathlib import PurePath
import gzip
import numpy as np
import dataclasses
from .file_buffers import (TwoLineFastaBuffer, FastQBuffer)
from .multiline_buffer import MultiLineFastaBuffer
from .bam import BamBuffer
from .delimited_buffers import (VCFBuffer, BedBuffer, GfaSequenceBuffer, get_bufferclass_for_datatype)
from .datatypes import GFFEntry, SAMEntry, ChromosomeSize
from .parser import NumpyFileReader, NpBufferedWriter, chunk_lines
from .chromosome_provider import FullChromosomeDictProvider, ChromosomeFileStreamProvider, LazyChromosomeDictProvider
from .indexed_fasta import IndexedFasta
from .npdataclassstream import NpDataclassStream
import logging

logger = logging.getLogger(__name__)
# from .gzip import gzip


class NpDataclassReader:
    def __init__(self, numpyfilereader):
        self._reader = numpyfilereader

    def __enter__(self):
        return self

    def __exit__(self):
        self._reader.close()

    def read(self):
        return self._reader.read().get_data()

    def read_chunk(self, chunk_size=5000000):
        chunk = self._reader.read_chunk(chunk_size)
        if chunk is None:
            # return an empty dataclass
            dataclass = self._reader._buffer_type.dataclass
            return dataclass(*[[] for field in dataclasses.fields(dataclass)])

        return chunk.get_data()

    def read_chunks(self, chunk_size=5000000):
        return NpDataclassStream((chunk.get_data() for chunk in self._reader.read_chunks(chunk_size)),
                                 dataclass=self._reader._buffer_type.dataclass)

    def __iter__(self):
        return self.read_chunks()


buffer_types = {
    ".vcf": VCFBuffer,
    ".bed": BedBuffer,
    ".fasta": MultiLineFastaBuffer,
    ".fa": MultiLineFastaBuffer,
    ".fastq": FastQBuffer,
    ".fq": FastQBuffer,
    ".gfa": GfaSequenceBuffer,
    ".gff": get_bufferclass_for_datatype(GFFEntry),
    ".gtf": get_bufferclass_for_datatype(GFFEntry),
    ".gff3": get_bufferclass_for_datatype(GFFEntry),
    ".sam": get_bufferclass_for_datatype(SAMEntry, comment="@"),
    ".bam": BamBuffer,
    ".sizes": get_bufferclass_for_datatype(ChromosomeSize)
}


def _get_buffered_file(
    filename, suffix, mode, is_gzip=False, buffer_type=None, **kwargs
):
    open_func = gzip.open if is_gzip else open
    if buffer_type is None:
        buffer_type = _get_buffer_type(suffix)
    if mode in ("w", "write", "wb"):
        return NpBufferedWriter(open_func(filename, "wb"), buffer_type)

    kwargs2 = {key: val for key, val in kwargs.items() if key in ["has_header"]}
    file_reader = NumpyFileReader(open_func(filename, "rb"), buffer_type, **kwargs2)
    if is_gzip:
        file_reader.set_prepend_mode()
    return NpDataclassReader(file_reader)


def _get_buffer_type(suffix):
    if suffix in buffer_types:
        return buffer_types[suffix]
    else:
        raise RuntimeError(f"File format {suffix} does not have a default buffer type. "
                           f"Specify buffer_type argument using get_bufferclass_for_datatype function or"
                           f"use one of {str(list(buffer_types.keys()))[1:-1]}")


def bnp_open(filename, mode=None, **kwargs):
    path = PurePath(filename)
    suffix = path.suffixes[-1]
    is_gzip = suffix in (".gz", ".bam")
    if suffix == ".gz":
        suffix = path.suffixes[-2]
    if suffix == ".fai":
        assert mode not in ("w", "write", "wb")
        return IndexedFasta(filename[:-4], **kwargs)
    return _get_buffered_file(filename, suffix, mode, is_gzip=is_gzip, **kwargs)


def count_entries(filename, buffer_type=None):
    logger.info(f"Counting entries in {filename}")
    path = PurePath(filename)
    suffix = path.suffixes[-1]
    is_gzip = suffix in (".gz", ".bam")
    if suffix == ".gz":
        suffix = path.suffixes[-2]
    open_func = gzip.open if is_gzip else open
    if buffer_type is None:
        buffer_type = _get_buffer_type(suffix)

    file_reader = NumpyFileReader(open_func(filename, "rb"), buffer_type)
    if is_gzip:
        file_reader.set_prepend_mode()
    chunk_counts = (chunk.count_entries() for chunk in file_reader.read_chunks())
    return sum(chunk_counts)
