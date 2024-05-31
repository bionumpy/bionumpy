from typing import List

import numpy as np
from ..encoded_array import OneToOneEncoding
from .exceptions import EncodingError


class AlphabetEncoding(OneToOneEncoding):
    """
    Encoding for an alphabet. The encoding is one-to-one and the alphabet is
    defined by the input string. The encoding is case-insensitive.
    """

    def __init__(self, alphabet: str):
        self._raw_alphabet = [c.upper() for c in alphabet]
        self._is_initialized = False
        self._alphabet_size = len(self._raw_alphabet)

    def _initialize(self, force=True):
        if self._is_initialized and not force:
            return
        alphabet = self._raw_alphabet
        self._alphabet = np.array([ord(c) for c in alphabet], dtype=np.uint8)
        lower_alphabet = (self._alphabet + ord("a")-ord("A"))
        self._alphabet = self._alphabet
        self._lookup = np.full(256, 255, dtype=np.uint8)
        self._lookup[self._alphabet] = np.arange(len(alphabet))
        self._lookup[lower_alphabet] = np.arange(len(alphabet))
        self._mask = np.zeros(256, dtype=bool)
        self._mask[self._alphabet] = True
        self._mask[lower_alphabet] = True
        self._is_initialized = True

    def _encode(self, byte_array):
        self._initialize()
        ret = self._lookup[byte_array]
        if np.any(ret >= self._alphabet_size):
            offset = np.flatnonzero(ret.ravel()==255)[0]
            try:
                tmp = [chr(c) for c in byte_array.ravel()[ret.ravel()==255]][:10]
            except:
                tmp = f'b{byte_array.ravel()[ret.ravel()==255]}'
            raise EncodingError(f"Error when encoding {''.join(chr(c) for c in byte_array.ravel()[0:100])} "
                                f"to {self.__class__.__name__}. Invalid character(s): "
                                f"{tmp}{[ord(c) for c in tmp]}", offset)
        return ret

    def _decode(self, encoded):
        self._initialize()
        array = np.asarray(encoded)
        return self._alphabet[array]

    @property
    def alphabet_size(self)->int:
        """
        Get the size of the alphabet

        Returns
        -------
        int
            The size of the alphabet

        """
        self._initialize()
        return self._alphabet.size

    def get_alphabet(self)-> List[str]:
        """
        Get the alphabet

        Returns
        -------
        list[str]
            The alphabet

        """
        self._initialize()
        return [chr(c) for c in self._alphabet]

    def get_labels(self) -> List[str]:
        return self.get_alphabet()

    def __str__(self):
        return f"""AlphabetEncoding('{"".join(self.get_alphabet())}')"""

    def __repr__(self):
        return f"""AlphabetEncoding('{"".join(self.get_alphabet())}')"""

    def __eq__(self, other):
        self._initialize()
        if not isinstance(other, AlphabetEncoding):
            return False
        other._initialize()
        if len(self._alphabet) != len(other._alphabet):
            return False
        return np.all(self._alphabet == other._alphabet)

    def __hash__(self):
        return hash(repr(self))


ACTGEncoding = AlphabetEncoding("ACTG")
ACGTEncoding = AlphabetEncoding("ACGT")
ACTGnEncoding = AlphabetEncoding("ACTGn")
ACGTnEncoding = AlphabetEncoding("ACGTn")
DigitEncoding = AlphabetEncoding("0123456789")
DNAEncoding = ACGTEncoding
ACUGEncoding = AlphabetEncoding("ACUG")
RNAENcoding = ACUGEncoding
AminoAcidEncoding = AlphabetEncoding('ACDEFGHIKLMNPQRSTVWY*')
BamEncoding = AlphabetEncoding("=ACMGRSVTWYHKDBN")
CigarOpEncoding = AlphabetEncoding("MIDNSHP=X")


class FlatAlphabetEncoding(AlphabetEncoding):
    def _encode(self, *args, **kwargs):
        return super()._encode(*args, **kwargs).ravel()


StrandEncoding = FlatAlphabetEncoding("+-.")


def get_alphabet_encodings():
    return [ACTGEncoding, ACGTEncoding, ACTGnEncoding, ACGTnEncoding, DigitEncoding,
            DNAEncoding, ACUGEncoding, RNAENcoding, AminoAcidEncoding,
            BamEncoding, CigarOpEncoding, StrandEncoding]
