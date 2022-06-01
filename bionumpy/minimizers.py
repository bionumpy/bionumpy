from .kmers import KmerEncoding
from .util import convolution
import numpy as np


@convolution
def get_minimizers(sequence, window_size, k):
    kmers = KmerEncoding(k)(np.lib.stride_tricks.sliding_window_view(sequence, k))
    minimizers = np.lib.stride_tricks.sliding_window_view(kmers, window_size-k+1).min(axis=-1)
    return minimizers
