import numpy as np
from .encodings.alphabet_encoding import AlphabetEncoding
from .encodings import BaseEncoding
from .kmers import KmerEncoding
from .sequences import as_sequence_array, Sequences


class DNAToProtein:
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    from_encoding = AlphabetEncoding("TCAG")
    to_encoding = BaseEncoding
    _lookup = np.array([ord(c) for c in amino_acids], dtype=np.uint8)

    def __getitem__(self, key):
        return self._lookup[key]


class WindowFunction:
    def windowed(self, sequences):
        sequences = as_sequence_array(sequences, encoding=self._encoding)
        assert np.all(sequences.shape.lengths % self.window_size == 0)
        tuples = sequences.ravel().reshape(-1, self.window_size)
        tuples.encoding = self._encoding
        new_data = self(tuples)
        return Sequences(new_data, sequences.shape.lengths // self.window_size, encoding=self._table.to_encoding)


class Translate(WindowFunction):
    def __init__(self, table=DNAToProtein()):
        self._table = table
        self._encoding = table.from_encoding

    @property
    def window_size(self):
        return 3# np.log(self._table._lookup.size, self._table.from_encoding.alphabet_size)

    def __call__(self, sequence):
        e = sequence.encoding
        sequence = sequence[..., ::-1]
        sequence.encoding = e
        kmer = KmerEncoding(self.window_size,  alphabet_encoding=self._encoding)(sequence)
        return self._table[kmer]
