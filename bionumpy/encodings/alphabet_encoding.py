import numpy as np
from .base_encoding import Encoding


class AlphabetEncoding(Encoding):

    def __init__(self, alphabet: str):
        alphabet = [c.lower() for c in alphabet]
        self._alphabet = np.array([ord(c) for c in alphabet], dtype=np.uint8)
        upper_alphabet = self._alphabet + ord("A")-ord("a")
        self._lookup = np.zeros(256, dtype=np.uint8)
        self._lookup[self._alphabet] = np.arange(len(alphabet))
        self._lookup[upper_alphabet] = np.arange(len(alphabet))
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

    def __eq__(self, other):
        if not isinstance(other, AlphabetEncoding):
            return False
        return np.all(self._alphabet == other._alphabet)


ACTGEncoding = AlphabetEncoding("ACTG")
ACTGnEncoding = AlphabetEncoding("ACTGn")
AminoAcidEncoding = AlphabetEncoding('ACDEFGHIKLMNPQRSTVWY')
