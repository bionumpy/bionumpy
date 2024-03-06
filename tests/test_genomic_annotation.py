from bionumpy.genomic_data import Genome


def test_genomic_annotation(data_path):
    g = Genome.from_file(data_path  / 'hg38.chrom.sizes')
    a = g.read_annotation(data_path / 'small_gff.gff3')
    print(a.genes.gene_id)
    print(a.transcripts.transcript_id)
    print(a.exons.exon_id)
