import bionumpy as bnp
from bionumpy.bnpdataclass import bnpdataclass


@bnpdataclass
class VCFEntry:
    chromosome: str
    position: int
    id: str
    ref_seq: str
    alt_seq: str
    quality: str
    filter: str
    info: str

class VCFBuffer(bnp.io.delimited_buffers.DelimitedBuffer):
    dataclass = VCFEntry

def unique_intersect(variants_filename: str, filename_b: str, chrom_sizes_file: str, out_filename: str):
    genome = bnp.Genome.from_file(chrom_sizes_file, filter_function=None)
    variants = bnp.open(variants_filename).read()
    genome_mask = genome.read_intervals(filename_b).get_mask()
    mask = genome_mask[genome.read_locations(variants_filename, has_numeric_chromosomes=True)]
    with bnp.open(out_filename, 'w') as f:
        f.write(variants[mask])


if __name__ == '__main__':
    import sys
    unique_intersect(*sys.argv[1:])
