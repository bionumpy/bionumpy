from bionumpy.encoded_array import as_encoded_array, EncodedArray
import bionumpy.encodings.alphabet_encoding as ae
from bionumpy.util.testing import assert_encoded_array_equal
import hypothesis.strategies as st
from hypothesis import given
import bionumpy as bnp
from .strategies import get_strategy_from_encoding
objs = (getattr(ae, name) for name in dir(ae) if not name.startswith("_"))
encodings = [obj for obj in objs if isinstance(obj, ae.AlphabetEncoding)]


@st.composite
def encoding_and_unencoded_data(draw):
    encoding = draw(st.sampled_from(encodings))
    return encoding, draw(get_strategy_from_encoding(encoding)())


@given(encoding_and_unencoded_data())
def test_encode_decode(data):
    encoding, unencoded = data
    lower_unencoded = unencoded.upper()
    unencoded = as_encoded_array(unencoded)
    encoded = encoding.encode(unencoded)
    decoded = EncodedArray(encoding.decode(encoded), bnp.encodings.BaseEncoding)
    assert_encoded_array_equal(decoded, as_encoded_array(lower_unencoded))
