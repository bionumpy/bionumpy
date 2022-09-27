from .encodings import AlphabetEncoding
import numpy as np
CigarOpEncoding = AlphabetEncoding("MIDNSHP=X")

class CigarEncoding:
    _mask = np.uint32(2**4-1)

    def decode(cls, encoded):
        lengths = encoded >> 4
        ops = encoded & cls._mask
        

def count_reference_length(cigars):
    symbol = cigars & np.uint32(2**4-1)
    lengths = (cigars >> 4)
    consuming = CigarOpEncoding.encode([ord(c) for c in"MDN=X"])
    mask = (symbol == consuming[0])
    for consuming_symbol in consuming[1:]:
        mask = mask | (symbol == consuming_symbol)
    return np.sum(mask*lengths, axis=-1).astype(int)

        
    

        
