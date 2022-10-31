import hypothesis.strategies as st
import numpy as np
from functools import partial

integers = partial(st.integers,
                   min_value=np.iinfo(np.int64).min+1,
                   max_value=np.iinfo(np.int64).max-1)

floats = partial(st.floats,
                 min_value=np.finfo(np.float64).min,
                 max_value=np.finfo(np.float64).max)

ascii_text = partial(st.text,
                     alphabet=st.characters(blacklist_characters="\t\n",
                                            min_codepoint=0, max_codepoint=127))


def get_strategy_from_encoding(encoding):
    whitelist = encoding.get_alphabet()
    whitelist = set(whitelist + [c.lower() for c in whitelist])
    whitelist = "".join(whitelist)
    return partial(st.text, alphabet=whitelist, min_size=1)
