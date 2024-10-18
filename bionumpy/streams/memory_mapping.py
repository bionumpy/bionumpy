import pickle
import warnings
from typing import Iterator, Callable

import numpy as np

from bionumpy import EncodedRaggedArray, EncodedArray
from bionumpy.encoded_array import Encoding


class MemMapEncodedRaggedArray:
    def __init__(self, basename):
        self._basename = basename

    @classmethod
    def load(self, basename: str) -> EncodedRaggedArray:
        '''
        Load a read only memory-mapped encoded ragged array.

        Parameters
        ----------
        basename: str
            The base name of the memory-mapped files.


        Returns
        -------
        EncodedRaggedArray
        '''
        data = np.memmap(f'{basename}_data.dat', dtype=np.uint8, mode='r')
        lengths = np.memmap(f'{basename}_lengths.dat', dtype=np.int32, mode='r')
        with open(f'{basename}_encoding.pkl', 'rb') as f:
            encoding = pickle.load(f)
        return EncodedRaggedArray(EncodedArray(data, encoding), lengths)

    @classmethod
    def create(cls, loader_creator: Callable[[], Iterator[EncodedRaggedArray]], basename) -> EncodedRaggedArray:
        '''
        Create a memory-mapped encoded ragged array.

        Parameters
        ----------
        loader_creator
        basename

        Returns
        -------
        EncodedRaggedArray
        '''
        warnings.warn(
            f"{cls.__name__} is in an experimental stage and may change in the future.",
            category=FutureWarning,
            stacklevel=2
        )
        total_sequence_length = 0
        n_sequences = 0
        encoding = None
        for sequences in loader_creator():
            n_sequences += len(sequences)
            total_sequence_length += sequences.size
            if encoding is None:
                encoding = sequences.encoding
            else:
                assert encoding == sequences.encoding, f'Expected {encoding} but got {sequences.encoding}'
        with open(f'{basename}_encoding.pkl', 'wb') as f:
            pickle.dump(encoding, f)

        data = np.memmap(f'{basename}_data.dat', dtype=np.uint8, mode='w+', shape=total_sequence_length)
        lengths = np.memmap(f'{basename}_lengths.dat', dtype=np.int32, mode='w+', shape=n_sequences)
        data_offset = 0
        length_offset = 0
        for sequences in loader_creator():
            data[data_offset:data_offset + sequences.size] = sequences.raw().ravel()
            data_offset += sequences.size

            lengths[length_offset:length_offset + len(sequences)] = sequences.lengths
            length_offset += len(sequences)
        data.flush()
        lengths.flush()
        return EncodedRaggedArray(EncodedArray(data, encoding), lengths)

