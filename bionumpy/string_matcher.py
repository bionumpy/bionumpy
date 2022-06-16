from bionumpy.rollable import RollableFunction
from bionumpy.sequences import as_sequence_array
import itertools
import numpy as np
import bionumpy as bnp
import re
from bionumpy.sequences import Sequence

class StringMatcher(RollableFunction):
    def __init__(self, matching_sequence, encoding):
        self._encoding = encoding
        self._matching_sequence_array = as_sequence_array(matching_sequence, encoding=encoding)

    @property
    def window_size(self):
        return len(self._matching_sequence_array)

    def __call__(self, sequence):
        return np.all(sequence == self._matching_sequence_array, axis=-1)


class FixedLenRegexMatcher(RollableFunction):
    def __init__(self, matching_regex, encoding):
        self._sub_matchers = construct_fixed_len_regex_matchers(matching_regex, encoding)

    @property
    def window_size(self):
        return self._sub_matchers[0].window_size

    def __call__(self, sequence):
        union_of_sub_matches = self._sub_matchers[0](sequence)
        for matcher in self._sub_matchers:
            union_of_sub_matches = np.logical_or(union_of_sub_matches, matcher(sequence))
        return union_of_sub_matches

        return np.all(sequence == self._matching_sequence_array, axis=-1)


class MaskedStringMatcher(RollableFunction):
    def __init__(self, matching_sequence_array, mask):
        #assert isinstance(matching_sequence_array, Sequence), type(matching_sequence_array)
        assert isinstance(mask, np.ndarray)
        assert matching_sequence_array.shape == mask.shape
        self._matching_sequence_array = matching_sequence_array
        self._mask = mask

    @property
    def window_size(self):
        return len(self._matching_sequence_array)

    def __call__(self, sequence):
        direct_match = ( sequence == self._matching_sequence_array )
        masked_or_match = np.logical_or(direct_match, self._mask)
        return np.all(masked_or_match, axis=-1)

def construct_fixed_len_regex_matchers(matching_regex : str, encoding):
    r = re.compile('\[[^\]]+\]')
    hit = r.search(matching_regex)
    if hit is None:
        return [construct_wildcard_matcher(matching_regex, encoding)]
    else:
        start, end = hit.span()
        pre, post = matching_regex[0: start], matching_regex[end:]
        return list(itertools.chain.from_iterable(
                [construct_fixed_len_regex_matchers(pre+symbol+post, encoding)
                for symbol in matching_regex[start+1 : end-1]]))


def construct_wildcard_matcher(matching_regex : str, encoding):
    mask = np.array( [symbol=='.' for symbol in matching_regex] )

    assert encoding == bnp.encodings.ACTGEncoding, "NotImplemented: Support for other encodings awaits a generic way to replace '.'with an arbitrary symbol supported by the encoding"
    base_seq = as_sequence_array( matching_regex.replace('.', 'A'), encoding=encoding )

    return MaskedStringMatcher(base_seq, mask)




