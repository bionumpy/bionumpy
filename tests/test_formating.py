from bionumpy.util.formating import table


def test_table():
    chromosomes = ["Chr10", "Chr20_test_test_test", "chr5"]
    sizes = [1, 100, 1000]
    t = table(zip(chromosomes, sizes))

    for value in chromosomes + sizes:
        assert str(value) in t