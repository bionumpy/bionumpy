from bionumpy.translate import Translate
from bionumpy import as_sequence_array


def test_translate():
    seqs = ["TTTTTC", "GGG"]
    prots = Translate().windowed(seqs)
    print(prots)
    assert prots == as_sequence_array(["FF", "G"])
