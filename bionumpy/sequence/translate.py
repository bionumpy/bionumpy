import numpy as np
from ..streams import streamable
from ..encodings import BaseEncoding, AlphabetEncoding
from ..sequence.kmers import KmerEncoder
from ..encoded_array import EncodedArray, as_encoded_array, EncodedRaggedArray
from ..bnpdataclass.bnpdataclassfunction import apply_to_npdataclass
from ..util.typing import EncodedArrayLike


class DNAToProtein:
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    from_encoding = AlphabetEncoding("TCAG") # AlphabetEncoding("TCAG")
    to_encoding = BaseEncoding
    _lookup = EncodedArray(np.array([ord(c) for c in amino_acids], dtype=np.uint8),
                           BaseEncoding)

    def __getitem__(self, key):
        return self._lookup[key.raw()]


class WindowFunction:
    def windowed(self, sequences):
        sequences = as_encoded_array(sequences, target_encoding=self._encoding)
        assert np.all(sequences.lengths % self.window_size == 0)
        tuples = sequences.ravel().reshape(-1, self.window_size)
        tuples.encoding = self._encoding
        new_data = self(tuples)
        return EncodedRaggedArray(EncodedArray(new_data, self._table.to_encoding),
                                  sequences.lengths // self.window_size)


class Translate(WindowFunction):
    def __init__(self, table=DNAToProtein()):
        self._table = table
        self._encoding = table.from_encoding

    @property
    def window_size(self):
        return 3# np.log(self._table._lookup.size, self._table.from_encoding.alphabet_size)

    def __call__(self, sequence: EncodedArrayLike):
        e = sequence.encoding
        sequence = sequence[..., ::-1]
        sequence.encoding = e
        kmer = KmerEncoder(self.window_size, alphabet_encoding=self._encoding)(sequence)
        return self._table[kmer]


@streamable()
@apply_to_npdataclass("sequence")
def translate_dna_to_protein(sequence):
    return Translate().windowed(sequence)
