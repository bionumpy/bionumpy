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


class DigitEncoding(Encoding):
    MIN_CODE = ord("0")

    @classmethod
    def encode(cls, bytes_array):
        if not isinstance(bytes_array, np.ndarray):
            bytes_array = bytes_array.raw()
        return bytes_array - cls.MIN_CODE

    @classmethod
    def decode(cls, digits):
        return digits + cls.MIN_CODE


class GenotypeEncoding(Encoding):
    @classmethod
    def encode(cls, bytes_array):
        assert bytes_array.shape[-1] == 3
        return (bytes_array[..., 0] == "1") + (
            bytes_array[..., 2] == "1"
        ).astype(np.int8)


class PhasedGenotypeEncoding:
    @classmethod
    def from_bytes(cls, bytes_array):
        assert bytes_array.shape[-1] == 3
        return 2*(bytes_array[..., 0] == ord("1")) + (
            bytes_array[..., 2] == ord("1")
        ).astype(np.int8)


class QualityEncoding(NumericEncoding):

    def encode(byte_array):
        assert np.all((byte_array >= ord("!")) & (byte_array < ord("!")+94)), repr(byte_array)
        res = byte_array - ord("!")
        return res

    def decode(quality):
        assert np.all(quality < 94), quality
        res = quality.astype(np.uint8) + ord("!")
        return res


def set_backend(lib):
    import sys

    from ..cupy_compatible.encodings.alphabet_encoding import CPAlphabetEncoding
    from ..cupy_compatible.encodings.alphabet_encoding import CPACTGEncoding
    from ..cupy_compatible.encodings.alphabet_encoding import CPAminoAcidEncoding
    
    sys.modules[__name__].AlphabetEncoding = CPAlphabetEncoding
    sys.modules[__name__].ACTGEncoding = CPACTGEncoding
    sys.modules[__name__].AminoAcidEncoding = CPAminoAcidEncoding
