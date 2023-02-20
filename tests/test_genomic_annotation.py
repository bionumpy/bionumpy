from bionumpy.genomic_data import Genome


def test_genomic_annotation():
    g = Genome.from_file('example_data/hg38.chrom.sizes')
    a = g.read_annotation('example_data/small_gff.gff3')
    print(a.genes.gene_id)
    print(a.transcripts.transcript_id)
    print(a.exons.exon_id)
