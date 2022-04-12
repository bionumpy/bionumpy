from pathlib import PurePath
import gzip
from .file_buffers import *
from .parser import *
from .delimited_buffers import *
from .chromosome_provider import *
from .indexed_fasta import IndexedFasta

buffer_types = {".vcf": VCFBuffer,
                ".bed": BedBuffer,
                ".fasta": TwoLineFastaBuffer,
                ".fa": TwoLineFastaBuffer,
                ".fastq": FastQBuffer,
                ".fq": FastQBuffer}

wrappers = {"chromosome_stream": ChromosomeStreamProvider,
            "dict": ChromosomeDictProvider,
            "stream": lambda x: x}

default_modes = {".vcf": "chromosome_stream",
                 ".bed": "dict"}

def get_buffered_file(filename, suffix, mode, is_gzip=False):
    open_func = gzip.open if is_gzip else open
    buffer_type = buffer_types[suffix]
    if mode in ("w", "write", "wb"):
        return NpBufferedWriter(open_func(filename, "wb"), buffer_type)

    buffers = NpBufferStream(open_func(filename, "rb"), buffer_type)
    data = (buf.get_data() for buf in buffers)
    if mode is None:
        mode = default_modes.get(suffix, "stream")
    return wrappers[mode](data)

def bnp_open(filename, mode=None, **kwargs):
    path = PurePath(filename)
    suffix = path.suffixes[-1]
    is_gzip = suffix == ".gz"
    if suffix == ".gz":
        suffix = path.suffixes[-2]
    if suffix in buffer_types:
        return get_buffered_file(filename, suffix, mode, is_gzip=is_gzip)
    if suffix == ".fai":
        assert mode not in ("w", "write", "wb")
        return IndexedFasta(filename[:-4], **kwargs)
