import numpy as np
from .encodings import ACTGTwoBitEncoding
from npstructures import RaggedArray
import logging


class KmerHash:
    def __init__(self, alphabet_size, k):
        self._powers = alphabet_size**np.arange(k)

    def hash(self, kmer):
        return kmer.dot(self._powers)


class TwoBitHash:
    def __init__(self, k=31, dtype=np.uint64):
        self._dtype = dtype
        self.k = k
        self._mask = dtype(4 ** k - 1)
        self._n_letters_in_dtype = 4 * dtype(0).nbytes
        self._shifts = dtype(2) * np.arange(self._n_letters_in_dtype, dtype=dtype)
        self._rev_shifts = self._shifts[::-1] + dtype(2)

    def get_kmers_with_buffer(self, sequence):
        res = sequence[:-1, None] >> self._shifts
        res |= sequence[1:, None] << self._rev_shifts
        return res

    def get_kmers(self, sequence, is_internal=False):
        assert sequence.dtype == self._dtype, "%s != %s" % (sequence.dtype, self._dtype)
        result = sequence[:, None] >> self._shifts
        result[:-1] |= sequence[1:, None] << self._rev_shifts
        return result

    def _has_buffer(self, sequences):
        last_end = sequences.shape.size  # intervals[1][-1]
        return (
            sequences._data.size * 4 - last_end >= self._n_letters_in_dtype - self.k + 1
        )

    def get_kmer_hashes(self, sequences):
        data = ACTGTwoBitEncoding.from_bytes(sequences._data)
        shape = sequences.shape
        func = (
            self.get_kmers_with_buffer
            if self._has_buffer(sequences)
            else self.get_kmers
        )
        kmers = func(data.view(self._dtype)).ravel()
        ra = RaggedArray(kmers, shape)

        try:
            ra = ra[:, : -(self.k - 1)]
        except IndexError:
            logging.error("Shape is %s" % (shape))
            raise
        return ra
        #return ACTGTwoBitEncoding.complement(ra._data) & self._dtype(4 ** self.k - 1)


class _KmerHash:
    CODES = np.zeros(256, dtype=np.uint64)
    CODES[[97, 65]] = 0
    CODES[[99, 67]] = 1
    CODES[[116, 84]] = 2
    CODES[[103, 71]] = 3

    def __init__(self, k=31):
        self.k = k
        self.POWER_ARRAY = np.power(4, np.arange(k, dtype=np.uint64), dtype=np.uint64)
        self.REV_POWER_ARRAY = np.power(
            4, np.arange(k, dtype=np.uint64)[::-1], dtype=np.uint64
        )

    def _to_text(self, seq):
        return "".join(seq[n] for n in seq)

    def _get_kmers(self, seq):
        return [
            self._to_text(seq[i : i + self.k]) for i in range(len(seq) - self.k + 1)
        ]

    def _get_codes(self, seq):
        return self.CODES[seq]

    def _join(self, kmers1, kmers2, mask):
        return np.concatenate((kmers1[mask], kmers2[mask]))

    def get_kmer_hashes(self, sequences):
        codes = self._get_codes(sequences.sequences)
        assert codes.dtype == np.uint64, codes.dtype
        kmers = np.convolve(codes, self.POWER_ARRAY, mode="valid")
        reverse_kmers = np.convolve((codes + 2) % 4, self.REV_POWER_ARRAY, mode="valid")
        mask = get_kmer_mask(
            (sequences.intervals_start, sequences.intervals_end),
            sequences.sequences.size,
            k=self.k,
        )
        return kmers, reverse_kmers, mask  # [mask], reverse_kmers[mask]
