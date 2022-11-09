from bionumpy.io.indexed_fasta import create_index
import bionumpy as bnp

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
