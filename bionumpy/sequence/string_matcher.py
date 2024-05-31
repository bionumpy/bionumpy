import logging
import itertools
import numpy as np
import re

from numpy.typing import ArrayLike

from .rollable import RollableFunction
from ..encoded_array import EncodedArray, as_encoded_array, Encoding
from ..util import as_strided
from npstructures import RaggedArray
from ..encodings import AlphabetEncoding
from ..util.typing import EncodedArrayLike, SingleEncodedArrayLike


def match_string(sequence: EncodedArrayLike, matching_sequence: SingleEncodedArrayLike) -> ArrayLike:
    """
    Matches a sequence aginst sequences and returns a boolean RaggedArray representing positions
    where the sequence matches.
    Parameters
    ----------
    sequence :
    matching_sequence :

    Returns
    -------
    ArrayLike
        A boolean RaggedArray representing positions where the sequence matches.

    Examples
    --------
    >>> import bionumpy as bnp
    >>> sequence = bnp.as_encoded_array(['ACGT', 'TACTAC'])
    >>> matching_sequence = bnp.as_encoded_array('AC', sequence.encoding)
    >>> bnp.match_string(sequence, matching_sequence)
    ragged_array([ True False False]
    [False  True False False  True])
    """
    sequence = as_encoded_array(sequence)
    enforced_encoding = sequence.encoding
    matching_sequence = as_encoded_array(matching_sequence, enforced_encoding)
    return StringMatcher(matching_sequence, enforced_encoding).rolling_window(sequence)


class StringMatcher(RollableFunction):
    def __init__(self, matching_sequence, encoding:Encoding):
        self._encoding = encoding
        self._matching_sequence_array = as_encoded_array(matching_sequence, target_encoding=encoding)

    @property
    def window_size(self):
        return len(self._matching_sequence_array)

    def __call__(self, sequence):
        return np.all(sequence == self._matching_sequence_array, axis=-1)

    
class RegexMatcher(RollableFunction):
    """
    Matches regexes of various lengths across a RaggedArray of sequences by constructing a list of FixedLenRegexMatcher objects from the original
    flexible length regex expression.

    It overrides the rolling_window function from the superclass to invoke FixedLenRegexMatcher objects across different window sizes for matcher
    objects.
    """

    def __init__(self, matching_regex, encoding):
        self._sub_matchers = construct_flexible_len_regex_matchers(matching_regex, encoding)

    def __call__(self, sequence: EncodedArray):
        raise NotImplementedError

    @property
    def window_size(self):
        return [sub_matcher.window_size for sub_matcher in self._sub_matchers]

    def rolling_window(self, _sequence: RaggedArray, window_size: int = None, mode="valid"):
        if hasattr(self, "_encoding") and self._encoding is not None:
            _sequence = as_encoded_array(_sequence, target_encoding=self._encoding)

        if mode == "valid":
            logging.warning("Mode is set to 'valid' in rolling_window(), but RegexMatcher uses only mode 'same'. Switching to 'same'...")

        shape, sequence = (_sequence.shape, _sequence.ravel())
        out = RaggedArray(np.zeros(sequence.shape, dtype=bool), shape)

        for index, sub_matcher in enumerate(self._sub_matchers):
            windows = as_strided(sequence, strides=sequence.strides + sequence.strides,
                                 shape=sequence.shape + (sub_matcher.window_size,), writeable=False)
            convoluted = sub_matcher(windows)

            if isinstance(_sequence, RaggedArray):
                out = np.logical_or(out, RaggedArray(convoluted, shape))
            elif isinstance(_sequence, (np.ndarray, EncodedArray)):
                out = np.logical_or(out, as_strided(convoluted, shape))

        return out


class FixedLenRegexMatcher(RollableFunction):
    def __init__(self, matching_regex, encoding):
        self._sub_matchers = construct_fixed_len_regex_matchers(matching_regex, encoding)
        self._encoding = encoding

    @property
    def window_size(self):
        return self._sub_matchers[0].window_size

    def __call__(self, sequence):
        union_of_sub_matches = self._sub_matchers[0](sequence)
        for matcher in self._sub_matchers:
            union_of_sub_matches = np.logical_or(union_of_sub_matches, matcher(sequence))
        return union_of_sub_matches


class MaskedStringMatcher(RollableFunction):
    def __init__(self, matching_sequence_array, mask):
        # assert isinstance(matching_sequence_array, Sequence), type(matching_sequence_array)
        assert isinstance(mask, np.ndarray)
        assert matching_sequence_array.shape == mask.shape
        self._matching_sequence_array = matching_sequence_array
        self._mask = mask

    @property
    def window_size(self):
        return len(self._matching_sequence_array)

    def __call__(self, sequence):
        assert sequence.shape[-1] == self.window_size, (sequence.shape, self._matching_sequence_array)
        direct_match = (sequence == self._matching_sequence_array)
        masked_or_match = np.logical_or(direct_match, self._mask)
        return np.all(masked_or_match, axis=-1)


def construct_fixed_len_regex_matchers(matching_regex: str, encoding):
    r = re.compile('\[[^\]]+\]')
    hit = r.search(matching_regex)
    if hit is None:
        return [construct_wildcard_matcher(matching_regex, encoding)]
    else:
        start, end = hit.span()
        pre, post = matching_regex[0: start], matching_regex[end:]
        return list(itertools.chain.from_iterable(
            [construct_fixed_len_regex_matchers(pre + symbol + post, encoding)
             for symbol in matching_regex[start + 1: end - 1]]))


def construct_flexible_len_regex_matchers(matching_regex: str, encoding):
    r = re.compile('(([A-Z]|\[[A-Z]+\])+)\.\{(\d*)\,(\d+)\}(.+)')
    hit = r.search(matching_regex)
    if hit is None:
        return construct_fixed_len_regex_matchers(matching_regex, encoding)
    else:

        min_gap = int(hit.group(3)) if hit.group(3) != '' else 0
        max_gap = int(hit.group(4))

        end_group_1 = hit.end(1)
        start_group_5 = hit.start(5)

        pre, post = matching_regex[0: end_group_1], matching_regex[start_group_5:]
        return list(itertools.chain.from_iterable(
            [construct_flexible_len_regex_matchers(pre + symbol + post, encoding)
             for symbol in [str("." * n) for n in range(min_gap, max_gap + 1)]]))


def construct_wildcard_matcher(matching_regex: str, encoding):
    assert isinstance(encoding, AlphabetEncoding)

    mask = np.array([symbol == '.' for symbol in matching_regex])
    #replacement = encoding.encoding.decode(0) if hasattr(encoding, "encoding") \
    #    else chr(encoding.decode(0))
    replacement = encoding._raw_alphabet[0]
    base_seq = as_encoded_array(matching_regex.replace('.', str(replacement)),
                                target_encoding=encoding)

    return MaskedStringMatcher(base_seq, mask)
