from pathlib import PurePath
import os
from .indexed_fasta import IndexedFasta, create_index
from .files import bnp_open
from .delimited_buffers import DelimitedBuffer
from .multiline_buffer import FastaIdx


class IndexBuffer(DelimitedBuffer):
    sep = "\t"
    dataclass = FastaIdx


def open_indexed(filename):
    path = PurePath(filename)
    suffix = path.suffixes[-1]
    #index_file_name = path+".fai2"
    index_file_name = path.with_suffix(path.suffix + ".fai")
    assert suffix in (".fa", ".fasta"), "Only fasta supported for indexed read"
    if not os.path.isfile(index_file_name):
        bnp_open(index_file_name, "w", buffer_type=IndexBuffer).write(create_index(path))
    return IndexedFasta(filename)
