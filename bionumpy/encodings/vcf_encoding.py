import itertools

import numpy as np

from ..encoded_array import EncodedArray, BaseEncoding
from bionumpy.encodings import Encoding


class _GenotypeRowEncoding(Encoding):
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

    def decode_lookup(self):
        _genotype_lookup = [
            sum([len(self._alphabet) ** (2 - i) * self._alphabet_lookup[ord(g)] for i, g in enumerate(genotype)])
            for genotype in self._genotypes
        ]
        _reverse_genotype_lookup = np.zeros((256, 4), dtype=np.uint8)
        _reverse_genotype_lookup[_genotype_lookup] = np.array([
            [ord(g) for g in genotype] + [ord("\t")] for genotype in self._genotypes
        ])

        return _reverse_genotype_lookup

    def encode(self, genotype_rows):
        data = self._preprocess_data_for_encoding(genotype_rows)
        n_rows = len(genotype_rows)
        encoded = \
            self._alphabet_size**2 * self._alphabet_lookup[data[:, 0].raw()] + \
            self._alphabet_size ** 1 * self._alphabet_lookup[data[:, 1].raw()] + \
            self._alphabet_size ** 0 * self._alphabet_lookup[data[:, 2].raw()]

        encoded = encoded.reshape(n_rows, len(encoded)//n_rows).astype(np.int8)
        return encoded

    def _preprocess_data_for_encoding(self, genotype_rows):
        if isinstance(genotype_rows, EncodedArray) and genotype_rows.ndim==3:
            return genotype_rows.reshape(-1, 3)
        # split the row of genotype data
        from ..io.strops import replace_inplace
        if isinstance(genotype_rows, list):
            assert len(genotype_rows) ==0
            genotype_rows = EncodedArray(np.zeros((0, 3)), BaseEncoding)
        data = genotype_rows.ravel()
        # hack because the row sometime ends with \n and sometimes with \t
        replace_inplace(data, "\n", "\t")
        indices = np.flatnonzero(data == "\t")
        indices = np.insert(indices, 0, -1)
        return data[indices[:-1, np.newaxis] + np.array([1, 2, 3])]

    def decode(self, genotype):
        if len(genotype.shape) == 0:
            return self.decode_lookup()[genotype]

        new_shape = genotype.shape[:-1] + (4*genotype.shape[-1],)
        if not isinstance(genotype, np.ndarray):
            genotype = genotype.raw()
        decoded = self.decode_lookup()[genotype].reshape(new_shape)
        # remove last tab
        return decoded[..., :-1]

    def to_string(self, e):
        if isinstance(e, np.ndarray):
            e = np.atleast_1d(e)

        if len(e.shape) == 2:
            return '\n'.join(self.to_string(c) for c in e)

        return ''.join(chr(c) for c in self.decode(e))

    def __repr__(self):
        return "GenotypeRowEncoding"


class GenotypeBuffer:
    _lookup = np.zeros(256, dtype=np.uint8)
    _lookup[[ord(c) for c in ('0', '1', '.')]] = np.array([0, 1, 0]) # Changed nan to 0 here since that's the uint8 version of nan

    def __init__(self, data, field_starts):
        self._data = data
        self._field_starts = field_starts
        self._genotype = None
        self._phasing_status = None

    def __getitem__(self, index):
         return self.__class__(self._data, self._field_starts[index])

    @property
    def genotype(self):
        if self._genotype is None:
            self._genotype = self._lookup[self._field_starts[..., None]+[0, 2]]
        return self._genotype

    @property
    def phasing_status(self):
        if self._phasing_status is None:
            self._phasing_status = self._lookup

        return


    def _preprocess_data_for_encoding(self, genotype_rows):
        # split the row of genotype data
        from ..io.strops import split, replace_inplace
        data = genotype_rows.ravel()
        # hack because the row sometime ends with \n and sometimes with \t
        replace_inplace(data, "\n", "\t")
        indices = np.flatnonzero(data == "\t")
        indices = np.insert(indices, 0, -1)
        return data[indices[:-1, np.newaxis] + np.array([1, 2, 3])]
        #data = split(data.ravel(), "\t")[:-1, 0:3]  # don't include last element which is empty
        #return data


class _PhasedGenotypeRowEncoding(_GenotypeRowEncoding):
    """Encoding that can be used when all records are phased
     and there is no missing data, i.e. every genotype is
     either 0|0, 0|1, 1|0 or 1|1.

     Does only support biallelic variants.
    """

    genotypes = ["0|0", "0|1", "1|0", "1|1"]

    def decode_lookup(self):
        return np.array([
            [ord(c) for c in genotype] + [ord("\t")]
            for genotype in self.genotypes], dtype=np.uint8)

    def encode(self, genotype_rows):
        if len(genotype_rows) == 0:
            return np.zeros((0, 1), dtype=np.uint8)
        data = self._preprocess_data_for_encoding(genotype_rows)
        n_rows = len(genotype_rows)
        encoded = (data[:, 0] == "1") * 2 + (data[:, 2] == "1")
        encoded = encoded.reshape(n_rows, -1).astype(np.int8)
        return encoded

    def __repr__(self):
        return "PhasedGenotypeRowEncoding"


class _PhasedHaplotypeRowEncoding(_GenotypeRowEncoding):
    """Encoding that encodes each haplotype (not two haplotypes together as _PhasdGenotypeRowEncoding"""
    _alleles = [str(i) for i in range(5)] + ["."]
    _alphabet = _alleles
    _reverse_alphabet_lookup = np.array([ord(c) for c in _alphabet], dtype=np.uint8)
    _alphabet_lookup = np.zeros(256, dtype=np.uint8)
    _alphabet_lookup[_reverse_alphabet_lookup] = np.arange(len(_reverse_alphabet_lookup))
    _alphabet_size = len(_alphabet)

    def encode(self, genotype_rows):
        data = self._preprocess_data_for_encoding(genotype_rows)
        n_rows = len(genotype_rows)
        if n_rows == 0:
            return np.zeros((0, 1), dtype=np.uint8)

        n_rows = len(genotype_rows)
        raw = data[:, 0].raw()
        first_haplotypes = self._alphabet_lookup[raw]
        second_haplotypes = self._alphabet_lookup[data[:, 2].raw()]
        out = np.zeros(len(first_haplotypes)*2, dtype=np.int8)
        out[::2] = first_haplotypes
        out[1::2] = second_haplotypes
        return out.reshape(n_rows, len(out)//n_rows)


PhasedGenotypeRowEncoding = _PhasedGenotypeRowEncoding()
PhasedHaplotypeRowEncoding = _PhasedHaplotypeRowEncoding()
GenotypeRowEncoding = _GenotypeRowEncoding()