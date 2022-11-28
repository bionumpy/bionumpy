import dataclasses
from itertools import chain
import numpy as np
from numpy.random import default_rng
from bionumpy.encoded_array import EncodedArray, as_encoded_array
from bionumpy.dna import reverse_compliment
from bionumpy.encodings import StrandEncoding

rng = default_rng()


@dataclasses.dataclass
class RNASeqSimulationSettings:
    transcript_counts: list = None
    fragment_size: int = 100
    sampling_rate: float = 0.9
    read_length: int = 75


def get_transcript_copies(sequences, sequence_counts):
    if sequence_counts is None:
        sequence_counts = [1] * len(sequences)
    indices = list(chain(*[[i] * count for i, count in enumerate(sequence_counts)]))
    return sequences[indices]


def fragment_transcript_copies(sequences, fragment_size):
    fragments = [sequence[i:i + fragment_size]
                 for sequence in sequences
                 for i in range(0, len(sequence) - fragment_size + 1, fragment_size)]
    return as_encoded_array(fragments, sequences.encoding)


def sample_transcript_fragments(sequences, sampling_rate):
    mask = np.random.choice(a=[True, False], size=len(sequences), p=[sampling_rate, 1 - sampling_rate])
    return sequences[mask]


def get_rnaseq_reads(fragments, read_length, strands=None):
    reverse_fragments = reverse_compliment(fragments)
    if strands is None:
        strands = EncodedArray(rng.choice([0, 1], replace=True, size=len(fragments)), StrandEncoding)
    reads = np.where(strands[:, np.newaxis] == "+", fragments[:, 0:read_length], reverse_fragments[:, 0:read_length])
    return reads


def simulate_rnaseq(reference_sequences, settings: RNASeqSimulationSettings):
    transcript_copies = get_transcript_copies(reference_sequences, settings.transcript_counts)
    fragmented_transcript_copies = fragment_transcript_copies(transcript_copies, settings.fragment_size)
    sampled_transcript_fragments = sample_transcript_fragments(fragmented_transcript_copies, settings.sampling_rate)
    rnaseq_reads = get_rnaseq_reads(sampled_transcript_fragments, settings.read_length)
    return rnaseq_reads
