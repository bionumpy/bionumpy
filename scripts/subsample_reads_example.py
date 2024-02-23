import dataclasses
from itertools import islice

import numpy as np
import bionumpy as bnp


@dataclasses.dataclass
class State:
    n_entries_left: int
    n_samples_left: int


def subsample_reads(input_filename: str, output_filename: str):
    bt = bnp.TwoLineFastaBuffer if input_filename.endswith(('fa', 'fasta', 'fa.gz', 'fasta.gz')) else None
    n_entries = bnp.count_entries(input_filename, buffer_type=bt)
    state = State(n_entries_left=n_entries,
                  n_samples_left=n_entries // 2)
    subsets = (sample_from_chunk(chunk, state) for chunk in bnp.open(input_filename, buffer_type=bt).read_chunks())
    with bnp.open(output_filename, "w", buffer_type=bt) as out_file:
        for subset in islice(subsets, 0, None):
            out_file.write(subset)


def sample_from_chunk(chunk: bnp.SequenceEntry, state: State) -> bnp.SequenceEntry:
    chunk_size = len(chunk)
    n_samples_in_this_chunk = np.random.hypergeometric(chunk_size, state.n_entries_left - chunk_size,
                                                       state.n_samples_left)
    rows = np.random.choice(np.arange(chunk_size), n_samples_in_this_chunk, replace=False)
    subset = chunk[np.sort(rows)]
    state.n_samples_left -= n_samples_in_this_chunk
    state.n_entries_left -= chunk_size
    return subset


def test():
    filename = 'example_data/big.fq.gz'
    out_filename = 'tmp.fq.gz'
    subsample_reads(filename, out_filename)
    assert bnp.count_entries(out_filename) == bnp.count_entries(filename) // 2


if __name__ == '__main__':
    import sys
    subsample_reads(sys.argv[1], sys.argv[2])
