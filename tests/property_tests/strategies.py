import hypothesis.strategies as st
import numpy as np
from functools import partial

integers = partial(st.integers,
                   min_value=np.iinfo(np.int64).min,
                   max_value=np.iinfo(np.int64).max)

def ascii_text():
    return st.text(alphabet=st.characters(min_codepoint=0, max_codepoint=127))
