from pathlib import PurePath
import gzip
from .file_buffers import (TwoLineFastaBuffer, FastQBuffer)
from .delimited_buffers import (VCFBuffer, BedBuffer, GfaSequenceBuffer, get_bufferclass_for_datatype)
from .datatypes import GFFEntry, SAMEntry
from .parser import NpBufferStream, NpBufferedWriter, chunk_lines
from .chromosome_provider import FullChromosomeDictProvider, ChromosomeFileStreamProvider
from .indexed_fasta import IndexedFasta

buffer_types = {
    ".vcf": VCFBuffer,
    ".bed": BedBuffer,
    ".fasta": TwoLineFastaBuffer,
    ".fa": TwoLineFastaBuffer,
    ".fastq": FastQBuffer,
    ".fq": FastQBuffer,
    ".gfa": GfaSequenceBuffer,
    ".gff": get_bufferclass_for_datatype(GFFEntry),
    ".gff3": get_bufferclass_for_datatype(GFFEntry),
    ".sam": get_bufferclass_for_datatype(SAMEntry)
}

wrappers = {
    "chromosome_stream": ChromosomeFileStreamProvider,
    "dict": FullChromosomeDictProvider,
    "stream": lambda x, y: x,
}

default_modes = {".vcf": "chromosome_stream", ".bed": "dict"}


def _get_buffered_file(
    filename, suffix, mode, is_gzip=False, buffer_type=None, **kwargs
):
    open_func = gzip.open if is_gzip else open
    if buffer_type is None:
        buffer_type = buffer_types[suffix]
    if mode in ("w", "write", "wb"):
        return NpBufferedWriter(open_func(filename, "wb"), buffer_type)

    kwargs2 = {key: val for key, val in kwargs.items() if key in ["chunk_size"]}
    buffers = NpBufferStream(open_func(filename, "rb"), buffer_type, **kwargs2)
    data = (buf.get_data() for buf in buffers)
    if "n_entries" in kwargs:
        data = chunk_lines(data, kwargs["n_entries"])
    if mode is None:
        mode = default_modes.get(suffix, "stream")
    return wrappers[mode](data, buffer_type.dataclass)


def bnp_open(filename, mode=None, **kwargs):
    path = PurePath(filename)
    suffix = path.suffixes[-1]
    is_gzip = suffix == ".gz"
    if suffix == ".gz":
        suffix = path.suffixes[-2]
    if suffix in buffer_types:
        return _get_buffered_file(filename, suffix, mode, is_gzip=is_gzip, **kwargs)
    if suffix == ".fai":
        assert mode not in ("w", "write", "wb")
        return IndexedFasta(filename[:-4], **kwargs)
