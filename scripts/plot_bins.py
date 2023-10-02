import bionumpy as bnp
import matplotlib

genome = bnp.genomic_data.BinnedGenome.from_file('example_data/sacCer3.fa', bin_size=10000)
genome.count_file('example_data/sacCer3_sample_variants.vcf.gz')
bnp.plot(genome).show()
