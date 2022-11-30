import bionumpy as bnp

rule translate_bionumpy:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/bionumpy/protein_sequences/{name}.fa"
    benchmark:
        "benchmarks/translate/bionumpy/{name}.txt"
    run:
        from bionumpy.translate import translate_dna_to_protein
        input_stream = bnp.open(input[0])
        output_stream = bnp.open(output[0], "w")
        output_stream.write(translate_dna_to_protein(input_stream))
        output_stream.close()

rule translate_biopython:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biopython/protein_sequences/{name}.fa"
    benchmark:
        "benchmarks/translate/biopython/{name}.txt"
    run:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        with open(output[0], 'w') as aa_fa:
            for dna_record in SeqIO.parse(input[0], 'fasta'):
                new_record = SeqRecord(
                    dna_record.seq.translate(),
                    id=dna_record.id)
                SeqIO.write(new_record, aa_fa, 'fasta')


rule translate_biostrings:
    input:
        "results/dna_sequences/{name}.fa"
    output:
        "results/biostrings/protein_sequences/{name}.fa"
    benchmark:
        "benchmarks/translate/biopython/{name}.txt"
    script:
        "scripts/reverse_complement_biostrings.R"
