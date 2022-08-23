from bionumpy.translate import Translate, DNAToProtein


def test_translate():
    seq = "TTTTTC"
    Translate(DNAToProtein())
