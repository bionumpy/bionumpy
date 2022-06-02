import numpy as np
import logging
from npstructures import RaggedArray
from .chromosome_map import ChromosomeMap

logger = logging.getLogger(__name__)


@ChromosomeMap()
def filter_on_intervals(entry, sorted_intervals):
    if len(sorted_intervals) == 0:
        mask = np.full(entry.position.shape, False)
    else:
        starts, ends = (sorted_intervals.start, sorted_intervals.end)
        idx = np.searchsorted(starts, entry.position, side="right") - 1
        idx = np.minimum(idx, starts.size - 1)
        mask = (entry.position >= starts[idx]) & (entry.position < ends[idx])
    return entry[mask]


@ChromosomeMap()
def get_snps(variants):
    snps = variants[variants.is_snp()]
    snps.ref_seq = snps.ref_seq.ravel()
    snps.alt_seq = snps.alt_seq.ravel()
    return snps


def convolution(func):
    def new_func(_sequence, window_size, *args, **kwargs):
        shape, sequence = (_sequence.shape, _sequence.ravel())
        convoluted = func(sequence, window_size, *args, **kwargs)
        if isinstance(_sequence, RaggedArray):
            out = RaggedArray(convoluted, shape)
        elif isinstance(_sequence, np.ndarray):
            out = np.lib.stride_tricks.as_strided(convoluted, shape)
        return out[..., : (-window_size + 1)]

    return new_func


def rolling_window_function(func):
    def new_func(_sequence, window_size, *args, **kwargs):
        shape, sequence = (_sequence.shape, _sequence.ravel())
        print(sequence, window_size)
        windows = np.lib.stride_tricks.sliding_window_view(sequence, window_size)
        convoluted = func(windows, window_size, *args, **kwargs)
        if isinstance(_sequence, RaggedArray):
            out = RaggedArray(convoluted, shape)
        elif isinstance(_sequence, np.ndarray):
            out = np.lib.stride_tricks.as_strided(convoluted, shape)
        return out[..., : (-window_size + 1)]

    return new_func


def pprint_one(sequence):
    return "".join(chr(c) for c in sequence)


def pprint(sequences):
    if isinstance(sequences, RaggedArray):
        return [pprint_one(s) for s in sequences]
    elif isinstance(sequences, np.ndarray):
        if len(sequences.shape) == 1:
            return pprint_one(sequences)
        return [pprint(s) for s in sequences]


def plot(obj):
    if not hasattr(obj, "__plot__"):
        logger.warning(f"{obj} has no __plot__ method")
