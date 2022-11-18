import itertools

import numpy as np
from bionumpy.encodings import Encoding


class GenotypeRowEncoding(Encoding):
    """
    Encoding for a row of genotype data in VCF-format.
    Supports phased and unphased genotypes and missing genotypes on the format ./.

    Does only support biallelic and triallelic variants.
    """

    # makes a flexible lookup by combining possible alleles
    # with the possible seperators. Can be extended by adding alleles or seperators
    _alleles = ["0", "1", "2", "."]
    _seperators = ["|", "/"]
    _alphabet = _alleles + _seperators
    _reverse_alphabet_lookup = np.array([ord(c) for c in _alphabet], dtype=np.uint8)
    _alphabet_lookup = np.zeros(256, dtype=np.uint8)
    _alphabet_lookup[_reverse_alphabet_lookup] = np.arange(len(_reverse_alphabet_lookup))
    _genotypes = list(itertools.product(_alleles, _seperators, _alleles))
    _alphabet_size = len(_alphabet)

    def get_labels(self):
        pass

    @classmethod
    def decode_lookup(cls):
        _genotype_lookup = [
            sum([len(cls._alphabet) ** (2 - i) * cls._alphabet_lookup[ord(g)] for i, g in enumerate(genotype)])
            for genotype in cls._genotypes
        ]
        _reverse_genotype_lookup = np.zeros((256, 4), dtype=np.uint8)
        _reverse_genotype_lookup[_genotype_lookup] = np.array([
            [ord(g) for g in genotype] + [ord("\t")] for genotype in cls._genotypes
        ])

        return _reverse_genotype_lookup

    @classmethod
    def encode(cls, genotype_rows):
        data = cls._preprocess_data_for_encoding(genotype_rows)
        n_rows = len(genotype_rows)
        encoded = \
            cls._alphabet_size**2 * cls._alphabet_lookup[data[:, 0].raw()] + \
            cls._alphabet_size ** 1 * cls._alphabet_lookup[data[:, 1].raw()] + \
            cls._alphabet_size ** 0 * cls._alphabet_lookup[data[:, 2].raw()]

        encoded = encoded.reshape(n_rows, len(encoded)//n_rows).astype(np.int8)
        return encoded

    @classmethod
    def _preprocess_data_for_encoding(cls, genotype_rows):
        # split the row of genotype data
        from ..io.strops import split, replace_inplace
        data = genotype_rows
        # hack because the row sometime ends with \n and sometimes with \t
        replace_inplace(data, "\n", "\t")
        data = split(data.ravel(), "\t")[:-1, 0:3]  # don't include last element which is empty
        return data

    @classmethod
    def decode(cls, genotype):
        if len(genotype.shape) == 0:
            return cls.decode_lookup()[genotype]

        new_shape = genotype.shape[:-1] + (4*genotype.shape[-1],)
        if not isinstance(genotype, np.ndarray):
            genotype = genotype.raw()
        decoded = cls.decode_lookup()[genotype].reshape(new_shape)# genotype.shape[0], genotype.shape[1]*4)
        # remove last tab
        return decoded[..., :-1]

    @classmethod
    def to_string(cls, e):
        if isinstance(e, np.ndarray):
            e = np.atleast_1d(e)

        if len(e.shape) == 2:
            return '\n'.join(cls.to_string(c) for c in e)

        return ''.join(chr(c) for c in cls.decode(e))


class PhasedGenotypeRowEncoding(GenotypeRowEncoding):
    """Encoding that can be used when all records are phased
     and there is no missing data, i.e. every genotype is
     either 0|0, 0|1, 1|0 or 1|1.

     Does only support biallelic variants.
    """

    genotypes = ["0|0", "0|1", "1|0", "1|1"]

    @classmethod
    def decode_lookup(cls):
        return np.array([
            [ord(c) for c in genotype] + [ord("\t")]
            for genotype in cls.genotypes], dtype=np.uint8)

    @classmethod
    def encode(cls, genotype_rows):
        data = cls._preprocess_data_for_encoding(genotype_rows)
        n_rows = len(genotype_rows)
        encoded = (data[:, 0] == "1") * 2 + (data[:, 2] == "1")
        encoded = encoded.reshape(n_rows, len(encoded)//n_rows).astype(np.int8)
        return encoded
