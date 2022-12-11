import biotite.sequence
import bionumpy as bnp


rule translate_bionumpy:
    input:
        "results/dna_sequences/{filename}.fa"
    output:
        "results/bionumpy/translate/{filename}.fa"
    benchmark:
        "benchmarks/translate/bionumpy/{filename}.txt"
    run:
        from bionumpy.sequence import translate_dna_to_protein
        input_stream = bnp.open(input[0]).read_chunks()
        output_stream = bnp.open(output[0], "w")
        output_stream.write(translate_dna_to_protein(input_stream))
        output_stream.close()

rule translate_biopython:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biopython/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/biopython/{name}.txt"
    script:
        "../scripts/biopython_translate.py"

rule translate_biostrings:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biostrings/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/biostrings/{name}.txt"
    script:
        "scripts/reverse_complement_biostrings.R"




rule translate_python:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/python/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/python/{name}.txt"
    script:
        "../scripts/python_translate.py"



rule translate_biotite:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biotite/translate/{name}.fa"
    benchmark:
        "benchmarks/translate/biotite/{name}.txt"

    run:
        import biotite.sequence.io.fasta as fasta
        fasta_file = fasta.FastaFile.read(input[0])
        out_fasta_file = fasta.FastaFile()

        for header, dna in fasta_file.items():
            dna = biotite.sequence.NucleotideSequence(dna)
            protein = dna.translate(complete=True)
            fasta.set_sequence(out_fasta_file, protein, header=header)

        out_fasta_file.write(output[0])
