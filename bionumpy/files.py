from pathlib import PurePath
import gzip
import numpy as np

from .file_buffers import (TwoLineFastaBuffer, FastQBuffer)
from .multiline_buffer import MultiLineFastaBuffer
from .delimited_buffers import (VCFBuffer, BedBuffer, GfaSequenceBuffer, get_bufferclass_for_datatype)
from .datatypes import GFFEntry, SAMEntry
from .parser import NumpyFileReader, NpBufferedWriter, chunk_lines
from .chromosome_provider import FullChromosomeDictProvider, ChromosomeFileStreamProvider, LazyChromosomeDictProvider
from .indexed_fasta import IndexedFasta
from .npdataclassstream import NpDataclassStream


class NpDataclassReader:
    def __init__(self, numpyfilereader):
        self._reader = numpyfilereader

    def __enter__(self):
        return self

    def __exit__(self):
        self._reader.close()

    def read(self):
        return self._reader.read().get_data()

    def read_chunks(self):
        return NpDataclassStream((chunk.get_data() for chunk in self._reader.read_chunks()), buffer_type=self._reader._buffer_type)

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
    ".sam": get_bufferclass_for_datatype(SAMEntry)
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
    is_gzip = suffix == ".gz"
    if suffix == ".gz":
        suffix = path.suffixes[-2]
    if suffix == ".fai":
        assert mode not in ("w", "write", "wb")
        return IndexedFasta(filename[:-4], **kwargs)
    return _get_buffered_file(filename, suffix, mode, is_gzip=is_gzip, **kwargs)
