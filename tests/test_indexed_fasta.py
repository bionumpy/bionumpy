from bionumpy.io.indexed_fasta import create_index
from bionumpy.datatypes import Interval
from bionumpy.util.testing import assert_encoded_array_equal
import bionumpy as bnp
import numpy as np


"""
0	300	3	80	81
1	600	310	80	81
2	900	921	80	81
3	1200	1836	80	81
"""


def test_fasta_index():
    index = create_index("example_data/small_genome.fa")


def test_dictlike():
    idx_fasta = bnp.open_indexed("example_data/small_genome.fa")
    assert list(idx_fasta.keys()) == ["0", "1", "2", "3"]
    assert "Indexed Fasta" in repr(idx_fasta)
    for key, val in idx_fasta.items():
        assert key in ["0", "1", "2", "3"]
        assert isinstance(val, bnp.EncodedArray)

    for val in idx_fasta.values():
        assert isinstance(val, bnp.EncodedArray)


def test_get_sequences():
    idx_fasta = bnp.open_indexed("example_data/small_genome.fa")
    intervals = Interval.from_entry_tuples([("1", 10, 20),
                                            ("2", 11, 50),
                                            ("1", 5, 10),
                                            ("3", 10, 110)])
    sequences = idx_fasta.get_interval_sequences(intervals)
    assert np.all(sequences.lengths == [10, 39, 5, 100])
    for interval, sequence in zip(intervals, sequences):
        print(interval)
        assert_encoded_array_equal(sequence, idx_fasta[interval.chromosome.to_string()][int(interval.start):int(interval.stop)])


