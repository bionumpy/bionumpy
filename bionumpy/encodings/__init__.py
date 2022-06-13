import numpy as np
from .base_encoding import BaseEncoding, Encoding
from .alphabet_encoding import AlphabetEncoding, ACTGEncoding, AminoAcidEncoding

__all__ = ["BaseEncoding", "Encoding",
           "AlphabetEncoding", "ACTGEncoding", "AminoAcidEncoding"]


class StrandEncoding(Encoding):
    MIN_CODE = ord("+")

    @classmethod
    def encode(cls, bytes_array):
        return (bytes_array & np.uint8(2)) >> np.uint8(1)

    @classmethod
    def decode(cls, strands):
        return 2 * strands + cls.MIN_CODE


class DigitEncoding(Encoding):
    MIN_CODE = ord("0")

    @classmethod
    def encode(cls, bytes_array):
        return bytes_array - cls.MIN_CODE

    @classmethod
    def decode(cls, digits):
        return digits + cls.MIN_CODE


class GenotypeEncoding(Encoding):
    @classmethod
    def encode(cls, bytes_array):
        assert bytes_array.shape[-1] == 3
        return (bytes_array[..., 0] == ord("1")) + (
            bytes_array[..., 2] == ord("1")
        ).astype("int")


class QualityEncoding(Encoding):
    def encode(byte_array):
        res = byte_array - ord("!")
        res.encoding = QualityEncoding
        return res

    def decode(quality):
        res = quality + ord("!")
        res.encoding = BaseEncoding
        return res
