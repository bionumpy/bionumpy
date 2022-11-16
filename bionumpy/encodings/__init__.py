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


class GenotypeRowEncoding(Encoding):
    lookup = np.array([
        [ord(chr) for chr in genotype]
        for genotype in ["0|0\t", "0|1\t", "1|0\t", "1|1\t"]], dtype=np.uint8)

    def get_labels(self):
        pass

    @classmethod
    def encode(cls, genotype_rows):
        from ..io.strops import split, replace_inplace
        data = genotype_rows
        n_rows = len(genotype_rows)
        # idea: Ravel, convert all, split on \t and reshape back to original shape
        # hack because the row sometime ends with \n and sometimes with \t
        replace_inplace(data, "\n", "\t")
        data = split(data.ravel(), "\t")[:-1, 0:3]  # don't include last element which is empty
        encoded = (data[:, 0] == "1") + 2 * (data[:, 2] == "1")
        encoded = encoded.reshape(n_rows, len(encoded)//n_rows).astype(np.int8)
        return encoded

    @classmethod
    def decode(cls, genotype):
        #assert genotype.shape[2] == 4, "Each decoded should be 4 bytes (e.g. 0|1\t)"
        new_shape = genotype.shape[:-1] + (4*genotype.shape[-1],)
        print(repr(genotype.data))
        decoded = cls.lookup[genotype.raw()].reshape(new_shape)# genotype.shape[0], genotype.shape[1]*4)
        # remove last tab
        return decoded[..., :-1]

    @classmethod
    def to_string(cls, e):
        return ''.join(chr(c) for c in cls.decode(np.atleast_1d(e)))


    """
    def encode(bytes_array):
        assert bytes_array.shape[-1] == 3
        return (bytes_array[..., 0] == "1") + (
            bytes_array[..., 2] == "1"
        ).astype(np.int8)
    """


"""
class PhasedGenotypeEncoding(GenotypeEncoding):
    #lookup = as_encoded_array(["0|0", "0|1", "1|0", "1|1"], BaseEncoding).to_numpy_array()
    lookup = np.array([
        [ord(chr) for chr in genotype]
        for genotype in ["0|0", "0|1", "1|0", "1|1"]], dtype=np.uint8)

    def encode(bytes_array):
        assert bytes_array.shape[-1] == 3
        return 2*(bytes_array[..., 0] == "1") + (
            bytes_array[..., 2] == "1"
        ).astype(np.int8)

    @classmethod
    def decode(cls, genotype):
        return cls.lookup[genotype]

"""



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

