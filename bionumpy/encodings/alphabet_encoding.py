import numpy as np
from ..encoded_array import OneToOneEncoding
from .exceptions import EncodingError


class AlphabetEncoding(OneToOneEncoding):
    def __init__(self, alphabet: str):
        self._raw_alphabet = [c.upper() for c in alphabet]
        self._is_initialized = False

    def _initialize(self):
        if self._is_initialized:
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
        if np.any(ret == 255):
            offset = np.flatnonzero(ret.ravel()==255)[0]
            tmp = [chr(c) for c in byte_array.ravel()[ret.ravel()==255]]
            raise EncodingError(f"Error when encoding {''.join(chr(c) for c in byte_array.ravel()[0:100])} "
                                f"to {self.__class__.__name__}. Invalid character(s): "
                                f"{tmp}", offset)
        return ret

    def _decode(self, encoded):
        self._initialize()
        return self._alphabet[np.asarray(encoded)]

    @property
    def alphabet_size(self):
        self._initialize()
        return self._alphabet.size

    def get_alphabet(self):
        self._initialize()
        return [chr(c) for c in self._alphabet]

    def get_labels(self):
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
AminoAcidEncoding = AlphabetEncoding('ACDEFGHIKLMNPQRSTVWY')
BamEncoding = AlphabetEncoding("=ACMGRSVTWYHKDBN")
CigarOpEncoding = AlphabetEncoding("MIDNSHP=X")
StrandEncoding = AlphabetEncoding("+-.")


def get_alphabet_encodings():
    return [ACTGEncoding, ACGTEncoding, ACTGnEncoding, ACGTnEncoding, DigitEncoding,
            DNAEncoding, ACUGEncoding, RNAENcoding, AminoAcidEncoding,
            BamEncoding, CigarOpEncoding, StrandEncoding]
