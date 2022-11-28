import numpy as np
from ..encoded_array import BaseEncoding, Encoding, NumericEncoding
from .alphabet_encoding import (AlphabetEncoding, DNAEncoding, RNAENcoding,
                                AminoAcidEncoding,
                                CigarOpEncoding, BamEncoding, StrandEncoding)


__all__ = ["BaseEncoding", "Encoding",
           "AlphabetEncoding", "ACTGEncoding", "AminoAcidEncoding"]# , "ACTGTwoBitEncoding"]


# class StrandEncoding(Encoding):
#     MIN_CODE = ord("+")
# 
#     @classmethod
#     def encode(cls, bytes_array):
#         return (bytes_array & np.uint8(2)) >> np.uint8(1)
# 
#     @classmethod
#     def decode(cls, strands):
#         return 2 * strands + cls.MIN_CODE


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
    #from ..cupy_compatible.encodings.alphabet_encoding import CPAlphabetEncoding
    #from ..cupy_compatible.encodings.alphabet_encoding import CPACTGEncoding
    #from ..cupy_compatible.encodings.alphabet_encoding import CPAminoAcidEncoding
    
    #sys.modules[__name__].AlphabetEncoding = CPAlphabetEncoding
    #sys.modules[__name__].ACTGEncoding = CPACTGEncoding
    #sys.modules[__name__].AminoAcidEncoding = CPAminoAcidEncoding

    from . import base_encoding
    base_encoding.np = lib

    from . import alphabet_encoding
    alphabet_encoding.np = lib

