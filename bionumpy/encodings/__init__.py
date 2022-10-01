import numpy as np
from ..sequences import create_sequence_array_from_already_encoded_data, ASCIIText
from .base_encoding import BaseEncoding, Encoding, NumericEncoding
from .alphabet_encoding import AlphabetEncoding, ACTGEncoding, AminoAcidEncoding
from ._legacy_encodings import ACTGTwoBitEncoding

__all__ = ["BaseEncoding", "Encoding",
           "AlphabetEncoding", "ACTGEncoding", "AminoAcidEncoding", "ACTGTwoBitEncoding"]


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
        return np.asarray(bytes_array) - cls.MIN_CODE

    @classmethod
    def decode(cls, digits):
        return digits + cls.MIN_CODE


class GenotypeEncoding(Encoding):
    @classmethod
    def encode(cls, bytes_array):
        assert bytes_array.shape[-1] == 3
        return (bytes_array[..., 0] == ord("1")) + (
            bytes_array[..., 2] == ord("1")
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
        byte_array = np.asarray(byte_array)
        assert np.all((byte_array >= ord("!")) & (byte_array < ord("!")+94)), repr(byte_array)
        res = byte_array - ord("!")
        return res

    def decode(quality):
        assert np.all(quality < 94)
        res = quality.astype(np.uint8) + ord("!")
        return create_sequence_array_from_already_encoded_data(res, ASCIIText)

