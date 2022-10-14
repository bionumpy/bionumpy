from .encodings.alphabet_encoding import CigarOpEncoding
from .encoded_array import as_encoded_array, EncodedArray
from npstructures import RaggedArray
import numpy as np


def split_cigar(cigars):
    if isinstance(cigars, RaggedArray):
        return map(lambda data: RaggedArray(data, cigars.shape), split_cigar(cigars.ravel()))

    symbol = EncodedArray((cigars & np.uint32(2**4-1)), CigarOpEncoding)
    lengths = (cigars >> 4)
    return symbol, lengths


def count_reference_length(symbol, lengths):
    consuming = as_encoded_array("MDN=X", CigarOpEncoding)
    mask = (symbol == consuming[0])
    for consuming_symbol in consuming[1:]:
        mask = mask | (symbol == consuming_symbol)

    return np.sum(mask*lengths, axis=-1).astype(int)

        
    

        
