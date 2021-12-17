import numpy as np
from .parser import get_mask_from_intervals
from .encodings import twobit_swap, ACTGTwoBitEncoding

def get_kmer_mask(intervals, size, k=31):
    starts, ends = intervals
    ends = np.maximum(starts, ends-k+1)
    return get_mask_from_intervals((starts, ends), size-k+1)


class TwoBitHash:
    def __init__(self, k=31, dtype=np.uint64):
        self._dtype = dtype
        self.k = k
        self._mask = dtype(4**k-1)
        self._n_letters_in_dtype = 4*dtype(0).nbytes
        self._shifts = dtype(2)*np.arange(self._n_letters_in_dtype, dtype=dtype)
        self._rev_shifts = self._shifts[::-1]+dtype(2)

    def _get_kmer_hashes(self, sequences):
        """Matching old interface"""
        last_end = sequences.intervals[1][-1]
        mask = get_kmer_mask(sequences.intervals, last_end, self.k)
        kmers = self.np_get_kmers(sequences.sequences.view(self._dtype))[:mask.size]
        reverse_hashes = sequences.encoding.complement(kmers) & self._dtype(4**self.k-1)
        forward_hashes = twobit_swap(kmers) >> ((self._n_letters_in_dtype-self.k)*2)
        return forward_hashes, reverse_hashes, mask

    def np_get_kmers_with_buffer(self, sequence):
        res = (sequence[:-1, None] >> self._shifts)
        res |= (sequence[1:, None] << self._rev_shifts)
        return res
    get_kmers_with_buffer = np_get_kmers_with_buffer

    def np_get_kmers(self, sequence, is_internal=False):
        assert sequence.dtype==self._dtype, sequence.dtype
        result = (sequence[:, None] >> self._shifts)
        result[:-1] |= (sequence[1:, None] << self._rev_shifts)
        if is_internal:
            return result
        result &= self._mask
        n_kmers = self._n_letters_in_dtype*sequence.size-self.k+1
        return result.ravel()[:n_kmers]

    def get_kmers(self, sequence):
        return self.np_get_kmers(sequence, is_internal=True)

    def _has_buffer(self, sequences):
        last_end = sequences.intervals[1][-1]
        return sequences.sequences.size*4 - last_end>=self._n_letters_in_dtype-self.k+1

    def get_new_kmer_hashes(self, sequences):
        last_end = sequences.intervals[1][-1]
        mask = get_kmer_mask(sequences.intervals, last_end, self.k)
        func = self.get_kmers_with_buffer if self._has_buffer(sequences) else self.get_kmers
        kmers = func(sequences.sequences.view(self._dtype))
        return kmers.ravel()[:mask.size][mask] & self._mask

    def get_kmer_hashes(self, sequences):
        return ACTGTwoBitEncoding.complement(self.get_new_kmer_hashes(sequences)) & self._dtype(4**self.k-1)

class KmerHash:
    CODES = np.zeros(256, dtype=np.uint64)
    CODES[[97, 65]] = 0
    CODES[[99, 67]] = 1
    CODES[[116, 84]] = 2
    CODES[[103, 71]] = 3

    def __init__(self, k=31):
        self.k=k
        self.POWER_ARRAY = np.power(4, np.arange(k, dtype=np.uint64), dtype=np.uint64)
        self.REV_POWER_ARRAY = np.power(4, np.arange(k, dtype=np.uint64)[::-1], dtype=np.uint64)

    def _to_text(self, seq):
        return "".join(letters[n] for n in seq)

    def _get_kmers(self, seq):
        return [self._to_text(seq[i:i+self.k]) for i in range(len(seq)-self.k+1)]

    def _get_codes(self, seq):
        return self.CODES[seq]

    def _join(self, kmers1, kmers2, mask):
        return np.concatenate((kmers1[mask], kmers2[mask]))

    def get_kmer_hashes(self, sequences):
        codes = self._get_codes(sequences.sequences)
        assert codes.dtype==np.uint64, codes.dtype
        kmers = np.convolve(codes, self.POWER_ARRAY, mode="valid")
        reverse_kmers = np.convolve((codes+2) % 4, self.REV_POWER_ARRAY, mode="valid")
        mask = get_kmer_mask((sequences.intervals_start, sequences.intervals_end), sequences.sequences.size, k=self.k)
        return kmers, reverse_kmers, mask # [mask], reverse_kmers[mask]
