import bionumpy as bnp


def analyze_across_chromosome(seq_fn, genes_fn):
    all_sequence = bnp.open(seq_fn).read_chunks(50)
    all_genes = bnp.open(genes_fn).read()
    genes_by_chr = bnp.groupby(all_genes, "chromosome")
    nn_counts = {nn:0 for nn in "ACTG"}
    for seq, (chr, genes) in zip(all_sequence, genes_by_chr):
        assert str(seq.name[0]) == chr, (seq, seq.name, chr)
        selected_seq = seq.sequence[0][genes]
        for nn in "ACTG":
            nn_counts[nn] += (selected_seq == nn).sum()
    gc_content = sum([nn_counts[nn] for nn in "GC"]) / sum(nn_counts.values())
    print(f"GC-content inside genes: {gc_content:.2f}", )

analyze_across_chromosome("example_data/gc_test_multichr.fa", "example_data/gc_bedtest_multichr.bed")
