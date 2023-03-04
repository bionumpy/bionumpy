import bionumpy as bnp
from bionumpy.genomic_data.geometry import Geometry


def test_pileup():
    chrom_sizes = bnp.open("example_data/hg38.chrom.sizes").read()
    geometry = Geometry.from_chrom_sizes(chrom_sizes)

    tf = bnp.open("example_data/ctcf.bed.gz").read()
    tf2 = bnp.open("example_data/znf263.bed.gz").read()
    
    
    reads = bnp.open("example_data/many_alignments.bam", buffer_type=bnp.io.bam.BamIntervalBuffer).read()
    print(set(reads.chromosome.tolist()))
    print(sum(c == 'chr1' for c in chrom_sizes.name.tolist()))
    pileup = geometry.get_pileup(reads)
    assert True
    # intervals_pileup = pileup.get_intervals(tf)


if __name__ == '__main__':
    # jaccard
    test_pileup()
