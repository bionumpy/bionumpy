import numpy as np

from bionumpy.genomic_data.binned_genome import BinnedGenome
from .genomic_fixtures import location_entries, genome_context, chrom_sizes


def test_count(location_entries, genome_context):
    bg = BinnedGenome(genome_context, bin_size=10)
    bg.count(location_entries)
    true = {
        'chr1': np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
        'chr2': np.array([0, 0, 1, 0, 0, 1, 0, 0, 0, 0])}
    for key in set(true.keys()) | set(bg.count_dict.keys()):
        np.testing.assert_equal(true[key], bg.count_dict[key])
