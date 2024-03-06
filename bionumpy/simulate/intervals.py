from typing import Dict

import numpy as np
from ..genomic_data import Genome
from ..datatypes import Interval
from .. import as_encoded_array


def simulate_fixed_size_uniform_intervals(genome: Genome, n_intervals: int, interval_size: int):
    chrom_sizes = genome.get_genome_context().chrom_sizes
    return simulate_intervals(chrom_sizes, interval_size, n_intervals)


def simulate_intervals(chrom_sizes: Dict[str, int], interval_size: int, n_intervals: int) -> Interval:
    names = as_encoded_array(list(chrom_sizes.keys()))
    sizes = np.array([i for i in chrom_sizes.values()])
    chromosome_probs = sizes / sizes.sum()
    simulated_chromosomes = np.random.choice(np.arange(len(sizes)), n_intervals, p=chromosome_probs)
    simulate_from_positions = np.zeros(n_intervals)
    simulate_to_positions = sizes[simulated_chromosomes] - interval_size
    start = np.random.randint(simulate_from_positions, simulate_to_positions)
    end = start + interval_size
    return Interval(names[simulated_chromosomes], start, end)
