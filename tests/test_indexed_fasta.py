from numpy.testing import assert_equal
from bionumpy.io.indexed_fasta import create_index
from bionumpy.datatypes import Interval
from bionumpy.util.testing import assert_encoded_array_equal
from bionumpy.encodings.string_encodings import StringEncoding
import bionumpy as bnp
import numpy as np


"""
0	300	3	80	81
1	600	310	80	81
2	900	921	80	81
3	1200	1836	80	81
"""


def test_fasta_index(data_path):
    index = create_index(data_path / "small_genome.fa")
    assert_equal(index.length, [300, 600, 900, 1200])


def test_dictlike(data_path):
    idx_fasta = bnp.open_indexed(data_path / "small_genome.fa")
    assert list(idx_fasta.keys()) == ["0", "1", "2", "3"]
    assert "Indexed Fasta" in repr(idx_fasta)
    for key, val in idx_fasta.items():
        assert key in ["0", "1", "2", "3"]
        assert isinstance(val, bnp.EncodedArray)

    for val in idx_fasta.values():
        assert isinstance(val, bnp.EncodedArray)


def test_get_sequences(data_path):
    idx_fasta = bnp.open_indexed(data_path / "small_genome.fa")
    _intervals = Interval.from_entry_tuples([("1", 10, 20),
                                            ("2", 11, 50),
                                            ("1", 5, 10),
                                             ("3", 10, 110),
                                             ("1", 80, 250)])
    intervals = bnp.bnpdataclass.replace(_intervals, chromosome=bnp.as_encoded_array(_intervals.chromosome, StringEncoding(['0', '1', '2', '3'])))
    sequences = idx_fasta.get_interval_sequences(intervals)
    assert np.all(sequences.lengths == [10, 39, 5, 100, 170])
    for interval, sequence in zip(intervals, sequences):
        assert sequence.to_string() == idx_fasta[interval.chromosome.to_string()][int(interval.start):int(interval.stop)].to_string()
