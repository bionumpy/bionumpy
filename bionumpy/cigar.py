from .encodings.alphabet_encoding import CigarOpArray
from .sequences import as_encoded_sequence_array
from npstructures import RaggedArray
import numpy as np


def split_cigar(cigars):
    if isinstance(cigars, RaggedArray):
        return map(lambda data: RaggedArray(data, cigars.shape), split_cigar(cigars.ravel()))

    symbol = (cigars & np.uint32(2**4-1)).view(CigarOpArray)
    lengths = (cigars >> 4)
    return symbol, lengths


def count_reference_length(symbol, lengths):
    consuming = as_encoded_sequence_array("MDN=X", CigarOpArray)
    mask = (symbol == consuming[0])
    for consuming_symbol in consuming[1:]:
        mask = mask | (symbol == consuming_symbol)

    return np.sum(mask*lengths, axis=-1).astype(int)

        
    

        
