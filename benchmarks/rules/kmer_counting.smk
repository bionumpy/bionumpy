rule jellyfish_count:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/jellyfish/kmer_counts/{filename}.jf"
    log:
        "{filename}.jf.log",
    params:
        kmer_length=8,
        size="1G",
        extra="--canonical",
    threads: 1
    wrapper:
        "v1.19.2-20-g6055e791/bio/jellyfish/count"

rule jellyfish_dump:
    input:
        "{prefix}.jf",
    output:
        "{prefix}.dump",
    log:
        "{prefix}.log",
    params:
        extra="-c -t",
    wrapper:
        "v1.19.2-20-g6055e791/bio/jellyfish/dump"


rule binoumpy_count:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/bionumpy/kmer_counts/{filename}.csv"
    run:
        from bionumpy.sequence import get_kmers, count_encoded
        from bionumpy.streams import streamable
        from bionumpy.io.dump_csv import dump_csv
        import bionumpy as bnp

        @streamable(sum)
        def count_kmers(sequence_entries):
            sequence = bnp.change_encoding(sequence_entries.sequence, bnp.DNAEncoding)
            kmers = get_kmers(sequence, k=8)
            return count_encoded(kmers, axis=None)

        stream = bnp.open(input[0]).read_chunks()
        output_stream = open(output[0], "wb")
        kmers = count_kmers(stream)
        output.write(bytes(dump_csv([(str, kmers.alphabet),
                                     (int, kmers.count)])))
                     

