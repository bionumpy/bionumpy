import gzip
import time
import bionumpy as bnp
import mgzip
import numpy as np
from isal import igzip

from bionumpy import MultiLineFastaBuffer
from bionumpy.io.npdataclassreader import NpDataclassReader
from bionumpy.io.parser import NumpyFileReader


def isal(file):
    n_chunks = 0
    with igzip.open(file, "rb") as f:
        while f.read(5000000):
            n_chunks += 1

    print("N chunks isal", n_chunks)

def np_filereader(file):
    with gzip.open(file, "rb") as f:
        return np.frombuffer(f.read(), dtype=np.uint8)


def np_frombuffer_chunks(file, chunk_size=5000000):
    n_chunks = 0
    with gzip.open(file, "rb") as f:
        while np.frombuffer(f.read(5_000_000), dtype=np.uint8).size:
            n_chunks += 1

    print("N chunks", n_chunks)


def np_filereader_read_chunks(file):
    with gzip.open(file, "rb") as f:
        reader = NumpyFileReader(f, buffer_type=bnp.io.FastQBuffer)
        reader.set_prepend_mode()
        chunks = reader.read_chunks(min_chunk_size=5000000)
        n_chunks = len([chunk.n_lines for chunk in chunks])

    print("N chunks np filereader", n_chunks)


def gzip_read_chunks(file):
    n = 0
    with gzip.open(file, "rb") as f:
        while f.read(5000000):
            pass

def read_python(file):
    with gzip.open(file, "rb") as f:
        return f.read()


def read_bionumpy(file):
    n_lines = 0
    with bnp.open(file) as f:
        for chunk in f.read_chunks(5000000):
            n_lines += len(chunk)
    print("N lines bionumpy", n_lines)


def read_mgzip(file):
    with mgzip.open(file, "rt", thread=8) as f:
        while f.read(5000000):
            pass


def read_npdataclass_with_isal(file):
    file_object = igzip.open(file, "rb")
    file_reader = NumpyFileReader(file_object, buffer_type=bnp.io.FastQBuffer)
    file_reader.set_prepend_mode()
    reader = NpDataclassReader(file_reader)
    n_lines = 0
    for chunk in reader.read_chunks():
        n_lines += len(chunk)
    print("N lines npdataclass isal", n_lines)

funcs = [isal, read_npdataclass_with_isal, gzip_read_chunks, read_bionumpy, np_filereader_read_chunks, np_filereader, read_python, read_mgzip]
#funcs = [np_filereader_read_chunks, np_frombuffer_chunks]
#funcs = [np_frombuffer_chunks]
#funcs = [np_filereader_read_chunks]


def benchmark(func, file):
    t0 = time.perf_counter()
    func(file)
    print("Time ", func.__name__, time.perf_counter() - t0)


if __name__ == "__main__":
    large_file = "../example_data/hg38.fa.gz"
    large_file = "../example_data/sacCer3.fa.gz"
    large_file = "../example_data/many_reads.fq.gz"

    for func in funcs:
        benchmark(func, large_file)

    t0 = time.perf_counter()
    n_lines = 0
    with bnp.open("../example_data/many_reads.fq") as f:
        for chunk in f:
            n_lines += len(chunk)

    print(time.perf_counter()-t0)


