import numpy as np
import dataclasses
from numpy.random import default_rng

from ..io.motifs import Motif
from ..sequence.position_weight_matrix import get_motif_scores
from ..datatypes import Interval, Bed6
from .. import streamable, EncodedArray, as_encoded_array
from ..encodings import AlphabetEncoding, StrandEncoding
from ..string_array import as_string_array

rng = default_rng()


@dataclasses.dataclass
class ChipSeqSimulationSettings:
    motif: Motif
    fragment_length: int = 200
    read_length: int = 100
    coverage: int = 10


def simulate_sequence(alphabet, length):
    numbers = rng.choice(np.arange(len(alphabet)), size=length)
    return EncodedArray(numbers, AlphabetEncoding(alphabet))


@streamable()
def simulate_chip_seq_fragments(reference_sequence, motif, n_fragments=1000, fragment_size=100):
    log_prob = get_motif_scores(reference_sequence, motif)# PWM.from_counts(motif))
    prob = np.exp(log_prob)
    prob /= prob.sum()
    points = rng.choice(np.arange(prob.size), size=n_fragments, replace=True, p=prob)
    left_extend = rng.poisson(fragment_size//2, size=points.size)
    right_extend = rng.poisson(fragment_size//2, size=points.size)
    start = np.maximum(points-left_extend, 0)
    stop = np.minimum(points+right_extend+1, log_prob.size)
    return Interval(["."]*len(start), start, stop)


@streamable()
def simulate_read_fragments(fragments: Interval, read_length: int):
    read_length = read_length
    strands = EncodedArray(rng.choice([0, 1], replace=True, size=len(fragments)), StrandEncoding)
    print(strands)
    starts = np.where(
        strands == "+",
        fragments.start,
        fragments.stop-read_length)
    stops = np.where(
        strands=="-",
        fragments.start+read_length,
        fragments.stop)
    starts = np.maximum(starts, fragments.start)
    stops = np.minimum(stops, fragments.stop)
    return Bed6(fragments.chromosome,
                starts,
                stops,
                ["."]*len(stops),
                [0]*len(stops),
                strands)


@streamable()
def simulate_chip_seq_reads(reference_sequence, settings, sequence_name=None):
    n_fragments = settings.coverage*len(reference_sequence) // settings.read_length
    fragments = simulate_chip_seq_fragments(reference_sequence, settings.motif,
                                            n_fragments, settings.fragment_length)
    reads = simulate_read_fragments(fragments, settings.read_length)
    if sequence_name is not None:
        print(repr(sequence_name))
        reads.chromosome = as_string_array([sequence_name]*len(reads))
    return reads
