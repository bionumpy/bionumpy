from typing import Dict

from bionumpy import as_encoded_array
from bionumpy.encodings.string_encodings import StringEncoding
from numpy.random import default_rng
import numpy as np
import npstructures as nps

from ..datatypes import SequenceEntry
from .. import EncodedRaggedArray, EncodedArray
from ..encodings.alphabet_encoding import AlphabetEncoding
from ..genomic_data.genomic_sequence import GenomicSequenceIndexedFasta, GenomicSequence
from ..datatypes import SequenceEntry, Interval, SequenceEntryWithQuality


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
 

def simulate_reads_from_genome(genome: GenomicSequence, length: int = 150, n_reads: int = 100,
                               chunk_size: int = 10000, sequence_name_prefix="", rng=default_rng(),
                               ignore_reads_with_n=False):
    """
    Simulates reads on a genome. Yields chunks of SequenceEntryWithQuality objects
    """
    if isinstance(rng, int):
        rng = default_rng(rng)

    chromosomes = genome.genome_context.chrom_sizes
    genome_size = sum([size for size in chromosomes.values()])

    for chromosome, chromosome_size in chromosomes.items():
        n_reads_on_chromosome = int(n_reads * chromosome_size / genome_size)
        n_simulated = 0
        while n_simulated < n_reads_on_chromosome:
            n_to_simulate = min(n_reads_on_chromosome - n_simulated, chunk_size)
            starts = rng.integers(0, chromosome_size-length, size=n_to_simulate)
            stops = starts + length

            chromosomes = as_encoded_array([chromosome] * n_to_simulate)
            intervals = Interval(chromosomes, starts, stops)
            sequences = genome.extract_intervals(intervals)

            names = as_encoded_array([f"{sequence_name_prefix}{i}" for i in range(n_simulated, n_simulated + n_to_simulate)])
            qualities = nps.RaggedArray(np.ones(sequences.size)*40, sequences.shape)
            sequence_entry = SequenceEntryWithQuality(
                names, sequences, qualities)

            if ignore_reads_with_n:
                n_mask = sequences == "N"
                n_mask = np.any(n_mask, axis=1)
                sequence_entry = sequence_entry[~n_mask]

            yield sequence_entry

            n_simulated += n_to_simulate
