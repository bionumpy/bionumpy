from typing import Union, List

from . import AlphabetEncoding
from ..encoded_array import Encoding
from ..encoded_array import EncodedArray, EncodedRaggedArray
from npstructures import RaggedArray
from ..util import is_subclass_or_instance
import numpy as np


class KmerEncoding(Encoding):
    """
    Encodes kmers of a given length using an alphabet encoding.
    """

    def __init__(self, alphabet_encoding: AlphabetEncoding, k: int):
        assert is_subclass_or_instance(alphabet_encoding, AlphabetEncoding), alphabet_encoding
        self._alphabet_encoding = alphabet_encoding
        self._k = k

    @property
    def k(self) -> int:
        return self._k

    def encode(self, data: Union[str, list, EncodedRaggedArray]) -> Union[EncodedArray, EncodedRaggedArray]:
        """
        Encodes a string, list of strings or EncodedRaggedArray into an EncodedArray or
        EncodedRaggedArray of hashed kmers.

        Parameters
        ----------
        data: Union[str, list, EncodedRaggedArray]

        Returns
        -------
        Union[EncodedArray, EncodedRaggedArray]

        """
        if isinstance(data, str):
            assert len(data) == self.k
            letters = self._alphabet_encoding.encode(data).raw()
            return EncodedArray(
                letters.dot(self._alphabet_encoding.alphabet_size ** np.arange(self._k)),
                self)
        if isinstance(data, (list, EncodedRaggedArray)):
            assert all(len(row) == self.k for row in data)
            letters = self._alphabet_encoding.encode(data).raw()
            if isinstance(letters, RaggedArray):
                letters = letters.to_numpy_array()
            return EncodedArray(
                letters.dot(self._alphabet_encoding.alphabet_size ** np.arange(self._k)),
                self)
        raise NotImplementedError

    def to_string(self, kmer: int) -> str:
        """
        Returns a human-readable string representation
        of an encoded kmer.
        """
        if np.asanyarray(kmer).ndim > 0:
            return ','.join(self.to_string(k) for k in kmer)
        if self._alphabet_encoding.alphabet_size == 4:
            tmp = (kmer >> (2 * np.arange(self._k))) & 3
        else:
            n = self._alphabet_encoding.alphabet_size
            tmp = (kmer // n ** np.arange(self._k)) % n

        chars = EncodedArray(tmp, self._alphabet_encoding)
        kmer = "".join(chr(int(b)) for b in chars.encoding.decode(chars).raw())
        return kmer

    def get_labels(self) -> List[str]:
        assert self._k <= 8, "Only supported for k <= 5"
        return [self.to_string(kmer) for kmer in range(self._alphabet_encoding.alphabet_size ** self._k)]

    def __str__(self):
        return f"{self._k}merEncoding({self._alphabet_encoding})"

    def __repr__(self):
        return f"KmerEncoding({self._alphabet_encoding}, {self._k})"

    def __eq__(self, other):
        if not isinstance(other, KmerEncoding):
            return False

        return self._k == other._k and self._alphabet_encoding == other._alphabet_encoding
