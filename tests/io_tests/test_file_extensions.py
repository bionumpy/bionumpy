import bionumpy as bnp

def test_fna(data_path):
    r = bnp.open(data_path / 'small.fna').read()
    assert len(r) == 3