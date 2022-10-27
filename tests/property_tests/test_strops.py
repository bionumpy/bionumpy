import pytest
from hypothesis import given, example, settings, Verbosity
import hypothesis.strategies as st
from bionumpy.io.strops import join, split, ints_to_strings
from numpy.testing import assert_array_equal
from bionumpy.testing import assert_encoded_array_equal, assert_encoded_raggedarray_equal
from bionumpy import as_encoded_array
from .strategies import ascii_text, integers

@pytest.mark.parametrize("sep", [",", "\t"])
@given(st.lists(ascii_text(), min_size=0))
def test_join(sep, strings):
    joined = join(as_encoded_array(strings), sep=sep)
    true = as_encoded_array(sep.join(strings))
    assert_encoded_array_equal(joined, true)

@given(ascii_text())
def test_split(sequence):
    seq = as_encoded_array(sequence)
    parts = split(seq, sep=",")
    assert_encoded_raggedarray_equal(
        parts,
        as_encoded_array(sequence.split(",")))

@given(st.lists(integers(), min_size=1))
def test_ints_to_strings(ints):
    strings = ints_to_strings(ints)
    assert_encoded_raggedarray_equal(
        strings,
        as_encoded_array([str(i) for i in ints], strings.encoding))
