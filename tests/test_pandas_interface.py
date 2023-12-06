import numpy as np
import pytest

from bionumpy.bnpdataclass.bnpdataclass import dynamic_concatenate
from bionumpy.util.testing import assert_bnpdataclass_equal
from .fixtures import data


@pytest.mark.xfail
@pytest.mark.parametrize('name', data.keys())
def test_pandas(name):
    d = data[name]
    df = d.topandas()
    d2 = d.__class__.from_data_frame(df)
    assert_bnpdataclass_equal(d, d2)
