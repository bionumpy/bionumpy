from ..encoded_array import BaseEncoding, Encoding, NumericEncoding
from .alphabet_encoding import (AlphabetEncoding, DNAEncoding, RNAENcoding,
                                AminoAcidEncoding, ACGTnEncoding,
                                CigarOpEncoding, BamEncoding, StrandEncoding)


__all__ = ["BaseEncoding", "Encoding",
           "AlphabetEncoding",  "AminoAcidEncoding"]


class DigitEncodingFactory(NumericEncoding):
    def __init__(self, min_code):
        self._min_code = ord(min_code)

    def _encode(self, bytes_array):
        return bytes_array - self._min_code

    def _decode(self, digits):
        return digits + self._min_code

    def __repr__(self):
        return "DigitEncoding(min_code=%d)" % self._min_code


DigitEncoding = DigitEncodingFactory("0")
QualityEncoding = DigitEncodingFactory("!")
CigarEncoding = DigitEncodingFactory(chr(0))


def set_backend(lib):

    from . import base_encoding
    base_encoding.np = lib

    from . import alphabet_encoding
    alphabet_encoding.np = lib

    from . import kmer_encodings
    kmer_encodings.np = lib

    from . import alphabet_encoding
    alphabet_encoding.np = lib

    # Have to reinitalize the encodings since these where initialized with numpy
    from .alphabet_encoding import DNAEncoding
    DNAEncoding._initialize(force=True)
