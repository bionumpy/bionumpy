from bionumpy.io.indexed_fasta import create_index


"""
0	300	3	80	81
1	600	310	80	81
2	900	921	80	81
3	1200	1836	80	81
"""


def test_fasta_index():
    index = create_index("example_data/small_genome.fa")
