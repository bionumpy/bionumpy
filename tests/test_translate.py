from bionumpy.sequence.translate import Translate
from bionumpy import as_encoded_array


def test_translate():
    seqs = ["TTTTTC", "GGG"]
    prots = Translate().windowed(seqs)
    print(prots)
    assert prots == as_encoded_array(["FF", "G"])
