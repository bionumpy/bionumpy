import numpy as np
from .base_encoding import Encoding
from ..sequences import ASCIIText, EncodedArray


class AlphabetEncoding(Encoding):
    def __init__(self, alphabet: str):
        alphabet = [c.lower() for c in alphabet]
        self._alphabet = np.array([ord(c) for c in alphabet], dtype=np.uint8)
        upper_alphabet = (self._alphabet + ord("A")-ord("a")).view(ASCIIText)
        self._alphabet = self._alphabet.view(ASCIIText)
        self._lookup = np.zeros(256, dtype=np.uint8)
        self._lookup[self._alphabet] = np.arange(len(alphabet))
        self._lookup[upper_alphabet] = np.arange(len(alphabet))
        # self._lookup.encoding = self
        self._mask = np.zeros(256, dtype=bool)
        self._mask[self._alphabet] = True
        self._mask[upper_alphabet] = True

    def encode(self, byte_array):
        return self._lookup[byte_array]

    def decode(self, encoded):
        return self._alphabet[np.asarray(encoded)]

    @property
    def alphabet_size(self):
        return self._alphabet.size

    def get_alphabet(self):
        return [chr(c) for c in self._alphabet]

    def __eq__(self, other):
        if not isinstance(other, AlphabetEncoding):
            return False
        return np.all(self._alphabet == other._alphabet)


def get_alphabet_array_class(alphabet):
    class AlphabetArray(EncodedArray):
        encoding = AlphabetEncoding(alphabet)

    AlphabetArray.__name__ = alphabet.upper()+"Array"
    AlphabetArray.__qualname__ = alphabet.upper()+"Array"
    return AlphabetArray


ACTGEncoding = AlphabetEncoding("ACTG")
ACTGnEncoding = AlphabetEncoding("ACTGn")
DNAEncoding = ACTGEncoding
ACUGEncoding = AlphabetEncoding("ACUG")
RNAENcoding = ACUGEncoding
AminoAcidEncoding = AlphabetEncoding('ACDEFGHIKLMNPQRSTVWY')
ACTGArray = get_alphabet_array_class("ACTG")
DNAArray = ACTGArray
RNAArray = get_alphabet_array_class("ACUG")
AminoAcidArray = get_alphabet_array_class('ACDEFGHIKLMNPQRSTVWY')
BamArray = get_alphabet_array_class("=ACMGRSVTWYHKDBN")
CigarOpArray = get_alphabet_array_class("MIDNSHP=X")
DigitArray = get_alphabet_array_class("0123456789")
