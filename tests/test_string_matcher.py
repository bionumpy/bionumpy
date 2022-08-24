import re

from npstructures import RaggedArray
from bionumpy.string_matcher import RegexMatcher
import bionumpy as bnp


def test_plain_string_matching():
    # Would probably reuse code for checking matches, but have different given, and different matcher class
    pass


def test_fixedlen_regex_matching():
    # TODO: These two should be given by hypothesis..
    sequences = ["ACGT", "AATGAT"]
    pattern = "[AG].[AT]"

    # to find (the fixed) pattern length
    pattern_len = len(re.sub("\[[^\]]+\]", ".", pattern))

    re_matcher = re.compile(pattern)

    re_matches = [[re_matcher.match(seq[offset:offset+pattern_len]) is not None
                   for offset in range(len(seq)-pattern_len+1)]
                  for seq in sequences]

    sequence_array = bnp.as_sequence_array(sequences, encoding=bnp.encodings.ACTGEncoding)
    matcher = RegexMatcher(pattern, encoding=bnp.encodings.ACTGEncoding, allow_flexible_lengths=False)
    matches = matcher.rolling_window(sequence_array)
    print(matches, re_matches)
    assert matches == RaggedArray(re_matches) #TODO: switch to correct assertion..


def test_flexible_len_regex_matching():
    sequences = ["ACGTTCG", "AATGATAA"]
    pattern = "AA.{,1}[CT]"

    re_matches = [[False for _ in range(5)],
                  [True, False, False, False, False, False]]

    sequence_array = bnp.as_sequence_array(sequences, encoding=bnp.encodings.ACTGEncoding)
    matcher = RegexMatcher(pattern, encoding=bnp.encodings.ACTGEncoding, allow_flexible_lengths=True)
    matches = matcher.rolling_window(sequence_array)
    print(matches, re_matches)
    assert matches == RaggedArray(re_matches)  # TODO: switch to correct assertion..
