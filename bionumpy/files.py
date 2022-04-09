from .file_buffers import *
from .parser import *
from .delimited_buffers import *
from .chromosome_provider import *
from .indexed_fasta import IndexedFasta

wrappers = {"chromosome_stream": ChromosomeStreamProvider,
            "dict": ChromosomeDictProvider,
            "stream": lambda x: x}

default_modes = {".vcf": "chromosome_stream",
                 ".bed": "dict"}

def get_buffered_file(filename, suffix, mode):
        buffers = BufferedNumpyParser.from_filename(filename)
        data = (buf.get_data() for buf in buffers)
        if mode is None:
            mode = default_modes.get(suffix, "stream")
        return wrappers[mode](data)

def bnp_open(filename, mode=None, **kwargs):
    path = PurePath(filename)
    suffix = path.suffixes[-1]
    if suffix == ".gz":
        suffix = path.suffixes[-2]
    if suffix in BufferedNumpyParser.buffer_types:
        return get_buffered_file(filename, suffix, mode)
    if suffix == ".fai":
        return IndexedFasta(filename[:-4], **kwargs)
