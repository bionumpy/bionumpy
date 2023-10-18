import dataclasses

import numpy as np
import bionumpy as bnp
from bionumpy.io.one_line_buffer import TwoLineFastaBuffer

snakemake = None


@dataclasses.dataclass
class State:
    n_entries_left: int
    n_samples_left: int


def subsample_reads(input_filename, output_filename):
    n_entries = bnp.count_entries(input_filename, buffer_type=TwoLineFastaBuffer)
    n_to_subsample = n_entries // 2
    n_entries_left = n_entries
    n_samples_left = n_to_subsample

    state = State(n_entries_left, n_samples_left)
    subsets = (sample_from_chunk(chunk, state) for chunk in bnp.open(input_filename, buffer_type=TwoLineFastaBuffer).read_chunks())

    with bnp.open(output_filename, "w", buffer_type=TwoLineFastaBuffer) as out_file:
        for subset in subsets:
            out_file.write(subset)
    assert state.n_samples_left==0
    assert state.n_entries_left==0


def sample_from_chunk(chunk, state):
    chunk_size = len(chunk)
    n_samples_in_this_chunk = get_number_of_samples(chunk_size, state.n_entries_left, state.n_samples_left)
    rows = np.random.choice(np.arange(chunk_size), n_samples_in_this_chunk, replace=False)
    rows.sort()
    subset = chunk[rows]
    state.n_samples_left -= n_samples_in_this_chunk
    state.n_entries_left -= chunk_size
    return subset


def get_number_of_samples(n_entries, n_entries_left, n_samples_left):
    if n_entries_left > n_entries:
        n_samples_in_this_chunk = np.random.hypergeometric(
            n_entries, n_entries_left - n_entries, n_samples_left)
    else:
        n_samples_in_this_chunk = n_samples_left
    return n_samples_in_this_chunk


def test_count_entries():
    filename = 'length150_nreads5000000.fa'
    bnp.count_entries(f'../benchmarks/results/dna_sequences/{filename}', buffer_type=TwoLineFastaBuffer)


def test():
    # filename = 'ENCFF689IPX.fq.gz'
    # filename = 'length150_nreads500000.fa'
    filename = 'length150_nreads5000000.fa'
    out = 'tmp.fa'
    input_filename = f'benchmarks/results/dna_sequences/{filename}'
    subsample_reads(input_filename, out)
    assert bnp.count_entries(out, buffer_type=TwoLineFastaBuffer) == bnp.count_entries(input_filename,
                                                                                       buffer_type=TwoLineFastaBuffer) // 2


if __name__ == '__main__':
    import sys

    subsample_reads(sys.argv[1], sys.argv[2])
