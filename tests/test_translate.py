from bionumpy.sequence.translate import Translate
from bionumpy.util.testing import assert_encoded_array_equal, assert_encoded_raggedarray_equal


def test_translate():
    seqs = ["TTTTTC", "GGG"]
    prots = Translate().windowed(seqs)
    assert_encoded_raggedarray_equal(prots, ["FF", "G"])
