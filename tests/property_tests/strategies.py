import hypothesis.strategies as st
import numpy as np
from functools import partial

integers = partial(st.integers,
                   min_value=np.iinfo(np.int64).min,
                   max_value=np.iinfo(np.int64).max)

floats = partial(st.floats,
                 min_value=np.finfo(np.float64).min,
                 max_value=np.finfo(np.float64).max)

ascii_text = partial(st.text,
                     alphabet=st.characters(blacklist_characters="\t\n",
                                            min_codepoint=0, max_codepoint=127))
