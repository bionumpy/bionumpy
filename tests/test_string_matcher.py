import re

from bionumpy.encodings import AlphabetEncoding
from npstructures import RaggedArray
from npstructures.testing import assert_raggedarray_equal

import bionumpy.encoded_array
from bionumpy.sequence.string_matcher import RegexMatcher, FixedLenRegexMatcher, StringMatcher
import bionumpy as bnp


def test_plain_string_matching():
    # Would probably reuse code for checking matches, but have different given, and different matcher class
    matcher = StringMatcher('V1-1', encoding=bnp.encodings.BaseEncoding)
    seq_array = bionumpy.encoded_array.as_encoded_array(['V1-1', 'V2-1', 'V1-1*2'])
    matches = matcher.rolling_window(seq_array)
    assert_raggedarray_equal(matches, RaggedArray([[True], [False], [True, False, False]]))


def test_fixedlen_regex_matching():
    # TODO: These two should be given by hypothesis..
    sequences = ["ACGT", "AATGAT"]
    pattern = "[AG].[AT]"

    # to find (the fixed) pattern length
    pattern_len = len(re.sub("\[[^\]]+\]", ".", pattern))

    re_matcher = re.compile(pattern)

    re_matches = [[re_matcher.match(seq[offset:offset + pattern_len]) is not None
                   for offset in range(len(seq) - pattern_len + 1)]
                  for seq in sequences]

    sequence_array = bionumpy.encoded_array.as_encoded_array(sequences, target_encoding=bnp.encodings.DNAEncoding)
    matcher = FixedLenRegexMatcher(pattern, encoding=bnp.encodings.DNAEncoding)
    matches = matcher.rolling_window(sequence_array)
    assert_raggedarray_equal(matches, RaggedArray(re_matches))  # TODO: switch to correct assertion..


def test_flexible_len_regex_matching():
    sequences = ["ACGTTCG", "AATGAAAC"]
    pattern = "AA.{,1}[CT]"

    re_matches = [[False for _ in range(7)],
                  [True, False, False, False, True, True, False, False]]

    sequence_array = bionumpy.encoded_array.as_encoded_array(sequences, target_encoding=bnp.encodings.DNAEncoding)
    matcher = RegexMatcher(pattern, encoding=bnp.encodings.DNAEncoding)
    matches = matcher.rolling_window(sequence_array)
    assert_raggedarray_equal(matches, RaggedArray(re_matches))  # TODO: switch to correct assertion..


def test_match_string():
    sequence_encoding = AlphabetEncoding("ACGT")
    sequences = bnp.as_encoded_array([
        "ACA", "TACTAC"
    ], sequence_encoding)

    correct = RaggedArray([
        [True, False],
        [False, True, False, False, True]
    ])
    pattern = "AC"
    matches = bnp.sequence.string_matcher.match_string(sequences, pattern)
    assert_raggedarray_equal(correct, matches)

