import bionumpy as bnp
from bionumpy.simulate.intervals import simulate_fixed_size_uniform_intervals


def _test_simulated_fixed_size_uniform_intervals():
    genome = bnp.Genome({"chr1": 1000, "chr15": 200})
    intervals = simulate_fixed_size_uniform_intervals(genome, 10, 20)
    assert len(intervals) == 10
    for interval in intervals:
        assert interval.chromosome.to_string() in ["chr1", "chr15"]
        assert interval.stop - interval.start == 20

