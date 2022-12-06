from typing import Dict
from numpy.random import default_rng
import numpy as np


from ..datatypes import SequenceEntry
from .. import EncodedRaggedArray, EncodedArray
from ..encodings.alphabet_encoding import AlphabetEncoding


def simulate_sequence(alphabet, length, rng=default_rng()):
    """Simulate a sequence from the given `alphabet` and of given `length`

    Parameters
    ----------
    alphabet : str
        The alphabet
    length : int
        The length of the sequence

    """
    numbers = rng.choice(np.arange(len(alphabet)), size=length)
    return EncodedArray(numbers, AlphabetEncoding(alphabet))


def simulate_sequences(alphabet: str, lengths: Dict[str, int], rng=default_rng()) -> SequenceEntry:
    """Simulate a set of sequences with the name and lengths of `lengths`

    Parameters
    ----------
    alphabet : str
        The alphabet to simulate from
    lengths : Dict[str, int]
        The names and lengths of the sequences

    Returns
    -------
    SequenceEntry
        SequenceEntry containing all the sequences

    """
    total_length = sum(lengths.values())
    flat_sequence = simulate_sequence(alphabet, total_length, rng=rng)
    names = list(lengths.keys())
    sequences = EncodedRaggedArray(flat_sequence, list(lengths.values()))
    se = SequenceEntry(names, sequences)
    return se
 
