import numpy as np
from .base_encoding import BaseEncoding, Encoding, NumericEncoding, CigarEncoding
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

    def encode(self, bytes_array):
        if not isinstance(bytes_array, np.ndarray):
            bytes_array = bytes_array.raw()
        return bytes_array - self._min_code

    def decode(self, digits):
        return digits + self._min_code


DigitEncoding = DigitEncodingFactory("0")
QualityEncoding = DigitEncodingFactory("!")



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

