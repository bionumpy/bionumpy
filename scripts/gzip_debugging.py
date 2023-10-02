import gzip
import time
import bionumpy as bnp
import mgzip
import numpy as np

from bionumpy import MultiLineFastaBuffer
from bionumpy.io.npdataclassreader import NpDataclassReader
from bionumpy.io.parser import NumpyFileReader


def np_filereader(file):
    with gzip.open(file, "rb") as f:
        return np.frombuffer(f.read(), dtype=np.uint8)
    
    
def np_filereader_read_chunks(file):
    with gzip.open(file, "rb") as f:
        return NumpyFileReader(f, buffer_type=MultiLineFastaBuffer).read_chunks()


def read_python(file):
    with gzip.open(file, "rb") as f:
        return f.read()


def read_bionumpy(file):
    with bnp.open(file) as f:
        return f.read()


def read_mgzip(file):
    with mgzip.open(file, "rt", thread=8) as f:
        return f.read()


funcs = [np_filereader_read_chunks, np_filereader, read_python, read_mgzip]


def benchmark(func, file):
    t0 = time.perf_counter()
    func(file)
    print("Time ", func.__name__, time.perf_counter() - t0)


if __name__ == "__main__":
    large_file = "../example_data/hg38.fa.gz"

    for func in funcs:
        benchmark(func, large_file)