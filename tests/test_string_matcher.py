import re

from npstructures import RaggedArray 
from bionumpy.string_matcher import FixedLenRegexMatcher
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
    matcher = FixedLenRegexMatcher(pattern, encoding=bnp.encodings.ACTGEncoding)
    matches = matcher.rolling_window(sequence_array)
    print(matches, re_matches)
    assert matches == RaggedArray(re_matches) #TODO: switch to correct assertion..

