from functools import reduce

import bionumpy as bnp


def unique_intersect(filename_a: str, filename_b: str, chrom_sizes_file: str, out_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_file, filter_function=None)
    intervals = genome.read_intervals(filename_a, buffer_type=bnp.Bed6Buffer)
    genome_mask = genome.read_intervals(filename_b).get_mask()
    mask = genome_mask[intervals].any(axis=-1)
    kept_intervals = intervals[mask]
    with bnp.open(out_filename, 'w') as f:
        f.write(kept_intervals.get_data())


def _test_profiling():
    name_a = 'ENCFF143HTO_mapped_reads_1m'
    name_b = 'ENCFF491EEI'
    unique_intersect(
                     f'../benchmarks/results/intervals/{name_a}.bed',
                     f'../benchmarks/results/bed_files/{name_b}.bed.gz',
                     '../example_data/hg38.chrom.sizes',
                     'tmp.bed')


# test_profiling()

if __name__ == '__main__':
    import sys
    unique_intersect(*sys.argv[1:])
