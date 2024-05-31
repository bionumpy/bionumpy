import bionumpy as bnp
import dataclasses


rule bionumpy_reverse_complement:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/bionumpy/reverse_complement/{filename}.fa"
    benchmark:
        repeat("benchmarks/reverse_complement/bionumpy/{filename}.txt", config["n_benchmark_repeats"])
    shell:
        "python3 ../scripts/reverse_compliment_example.py {input} {output}"
    # run:
    #     sequence_entries = bnp.open(input[0], buffer_type=bnp.TwoLineFastaBuffer).read_chunks()
    #     reversed_entries = bnp.sequence.get_reverse_complement(sequence_entries)
    #     with bnp.open(output[0], "w", buffer_type=bnp.TwoLineFastaBuffer) as outfile:
    #         outfile.write(reversed_entries)


rule seqtk_reverse_complement:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/seqtk/reverse_complement/{filename}.fa"
    log:
        "results/seqtk/reverse_complement/{filename}.log"
    params:
        extra="-r",
    benchmark:
        repeat("benchmarks/reverse_complement/seqtk/{filename}.txt", config["n_benchmark_repeats"])
    wrapper:
        "v1.19.2/bio/seqtk/seq"


rule biopython_reverse_complement:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/biopython/reverse_complement/{filename}.fa"
    benchmark:
        repeat("benchmarks/reverse_complement/biopython/{filename}.txt", config["n_benchmark_repeats"])
    conda:
        "../envs/biopython.yml"
    script: "../scripts/biopython_reverse_complement.py"


rule python_reverse_complement:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/python/reverse_complement/{filename}.fa"
    benchmark:
        repeat("benchmarks/reverse_complement/python/{filename}.txt", config["n_benchmark_repeats"])
    run:
        def reverse_complement(sequence):
            mapping = {"A": "T", "T": "A", "C": "G", "G": "C",
                       "a": "t", "t": "a", "c": "g", "g": "c"}
            reverse = sequence[::-1]
            complement = "".join(mapping[b] for b in reverse)
            return complement

        with open(input[0]) as infile:
            with open(output[0], "w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        outfile.write(line)
                    else:
                        outfile.write(reverse_complement(line.strip()) + "\n")


rule biotite_reverse_complement:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biotite/reverse_complement/{name}.fa"
    benchmark:
        repeat("benchmarks/reverse_complement/biotite/{name}.txt", config["n_benchmark_repeats"])
    conda:
        "../envs/biotite.yml"
    script:
        "../scripts/biotite_reverse_complement.py"
