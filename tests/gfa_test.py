from bionumpy.io.gfa import GfaPathBuffer
import pytest
import bionumpy as bnp


@pytest.mark.skip("Hard coded paths")
def test():
    f = bnp.open("/home/knut/Data/small_path_gfa.gfa", buffer_type=GfaPathBuffer)
    data = f.read()
    print(data.node_ids)
    assert False
