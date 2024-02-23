import bionumpy as bnp


def gc_genomic(sequence_filename: str, genes_filename: str) -> float:
    genome = bnp.Genome.from_file(sequence_filename)
    gene_mask = genome.read_intervals(genes_filename).get_mask()
    sequence = genome.read_sequence()
    gene_sequences = sequence[gene_mask]
    counts = bnp.count_encoded(gene_sequences, axis=None)
    gc_content = (counts["G"] + counts["C"]) / (gene_sequences.size - counts["N"])
    return gc_content


def test():
    assert gc_genomic("example_data/gc_test_multichr.fa", "example_data/gc_bedtest_multichr.bed") == 1 / 3


if __name__ == '__main__':
    import sys

    print(f'GC content {gc_genomic(*sys.argv[1:])}')
