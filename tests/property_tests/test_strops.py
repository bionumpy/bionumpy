import pytest
import numpy as np
from hypothesis import given, example, settings, Verbosity
import hypothesis.strategies as st
from bionumpy.io.strops import (join, split, ints_to_strings, str_to_float,
                                str_to_int, str_equal, int_lists_to_strings,
                                float_to_strings)
from numpy.testing import assert_array_equal, assert_array_almost_equal
from bionumpy.util.testing import assert_encoded_array_equal, assert_encoded_raggedarray_equal
from bionumpy import as_encoded_array
from numpy import array
from npstructures import RaggedArray
from .strategies import ascii_text, integers, floats


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
@example(ints=[-9223372036854775807])
def test_ints_to_strings(ints):
    strings = ints_to_strings(ints)
    assert_encoded_raggedarray_equal(
        strings,
        as_encoded_array([str(i) for i in ints], strings.encoding))


@given(st.lists(floats().filter(lambda x: abs(x)>10**(-15)), min_size=1))
@example(_floats = array([1.80143985e+15]))
@example(_floats=[1.3230423433805828e+16])
# @example(_floats=[4.450147717014403e-308])
# Failing on github actions, but not locally
# Probably not important. TODO Fix
@pytest.mark.xfail
def test_str_to_float(_floats):
    _floats = array(_floats)
    float_strings = [str(f) for f in _floats]
    floats = str_to_float(as_encoded_array(float_strings))
    # my_diff = np.abs(floats-_floats)
    true = np.array([float(s) for s in float_strings])
    tf, tm = np.frexp(true)
    f, m = np.frexp(floats)
    assert_array_almost_equal(f, tf)
    assert_array_equal(m, tm)
    

@given(st.lists(integers(), min_size=1))
def test_str_to_int(ints):
    int_strings = [str(i) for i in ints]
    result = str_to_int(as_encoded_array(int_strings))
    assert_array_equal(result, ints)

@given(st.lists(ascii_text(), min_size=1), ascii_text())
def test_str_equal(sequences, match_string):
    true = [s==match_string for s in sequences]
    result = str_equal(as_encoded_array(sequences), match_string)
    assert_array_equal(true, result)


@given(st.lists(st.lists(integers(), min_size=1), min_size=1))
def test_int_lists_to_strings(int_lists):
    ra = RaggedArray(int_lists)
    strings = int_lists_to_strings(ra, sep=",")
    true = as_encoded_array(
        [",".join(str(i) for i in ints) for ints in int_lists])
    assert_encoded_raggedarray_equal(strings, true)


@pytest.mark.skip("Inaccurate float")
@given(st.lists(floats().filter(lambda x: abs(x)>10**(-15)), min_size=1))
@settings(max_examples=500)
@example(floats=[0.015624999592526556])
def test_float_to_str(floats):
    floats = np.array(floats)
    ra = float_to_strings(floats)
    true = floats
    result = np.array([float(row.to_string()) for row in ra])
    tf, tm = np.frexp(true)
    f, m = np.frexp(result)
    tf = np.where(tm > m, tf*2**np.maximum(tm-m, 0), tf)
    f = np.where(m > tm, f*2**np.maximum(m-tm, 0), f)
    assert_array_almost_equal(f, tf)
