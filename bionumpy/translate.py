from npstructures import RaggedArray
import numpy as np
from .rollable import RollableFunction
from .encodings.alphabet_encoding import AlphabetEncoding
from .kmers import KmerEncoding
from .sequences import as_sequence_array

class DNAToProtein:
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    from_encoding = AlphabetEncoding("TCAG")
    _lookup = np.array([ord(c) for c in amino_acids], dtype=np.uint8)

    def __getitem__(self, key):
        return self._lookup[key]
    

class WindowFunction:
    def windowed(self, sequences):
        assert np.all(sequences.shape.lengths % self.window_size == 0)
        sequences = as_sequence_array(sequences)
        tuples = seqeunces.ravel().reshape(-1, self.window_size)
        new_data = self(tuples)
        return RaggedArray(new_data, sequences.shape.lengths // self.window_size)

class Translate(WindowFunction):
    def __init__(self, table=DNAToProtein()):
        self._table = table
        self._encoding=table.from_encoding

    @property
    def window_size(self):
        return np.log(self._table._lookup.size, self._table.from_encoding.alphabet_size)

    def __call__(self, sequence):
        kmer = KmerEncoding(self._table, alphabet_encoding=self._encoding)(sequence)
        return self._table[kmer]
        return protein_codes
